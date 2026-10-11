import copy
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest
import numpy as np
from reducao import presolve,check_certificates,fractional_bound,core,reconstruct,core_qubo,solve_core_qaoa,exact
from recursivo import eliminate,rqaoa,qiro,cluster_representatives


def enumerate_knapsack(c,v,B,fixed=None):
    fixed=fixed or {}
    feasible=[(list(b),sum((vi*bi for vi,bi in zip(v,b)),F(0))) for b in product((0,1),repeat=len(c))
              if sum(ci*bi for ci,bi in zip(c,b))<=B and all(b[i]==value for i,value in fixed.items())]
    return feasible

class TestReducao(unittest.TestCase):
    def test_certificados_preservam_todos_otimos_inclusive_empates(self):
        rng=np.random.default_rng(82)
        for _ in range(120):
            c=[int(x) for x in rng.integers(1,9,6)]
            v=[F(int(x),17) for x in rng.integers(0,15,6)];B=int(rng.integers(0,25))
            feasible=enumerate_knapsack(c,v,B);opt=max(value for b,value in feasible)
            result=presolve(c,v,B)
            self.assertTrue(check_certificates(c,v,B,result))
            for b,value in feasible:
                if value==opt:self.assertTrue(all(b[i]==bit for i,bit in result['fixed'].items()))
            _,dp=exact(c,v,B);self.assertEqual(dp,opt)

    def test_limite_fracionario_realmente_superior(self):
        rng=np.random.default_rng(13)
        for _ in range(60):
            c=[int(x) for x in rng.integers(1,7,5)];v=[F(int(x),19) for x in rng.integers(0,20,5)]
            B=int(rng.integers(0,15));fixed={0:int(rng.integers(0,2)),3:int(rng.integers(0,2))}
            feasible=enumerate_knapsack(c,v,B,fixed);ub=fractional_bound(c,v,B,fixed)
            if feasible:self.assertGreaterEqual(ub,max(value for b,value in feasible))
            else:self.assertIsNone(ub)

    def test_certificado_adulterado_rejeitado(self):
        c=[2,3,10];v=[F(2),F(1),F(3)];r=presolve(c,v,2)
        r['certificados'][0]['limite_inferior']['numerador']+=1
        self.assertFalse(check_certificates(c,v,2,r))

    def test_dominancia_simples_nao_elimina_controle_co_selecionavel(self):
        c=[1,2];v=[F(2),F(1)];r=presolve(c,v,3)
        self.assertEqual(r['fixed'],{0:1,1:1})

    def test_mdc_reconstrucao_offset_qubo_e_ising(self):
        c=[6,4,8,10];v=[F(7,10),F(3,10),F(9,10),F(1)];B=19;fixed={0:1}
        m=core(c,v,B,fixed);self.assertEqual(m['gcd'],2)
        for reduced in product((0,1),repeat=len(m['free'])):
            original=reconstruct(reduced,m,len(c))
            feasible=sum(ci*bi for ci,bi in zip(c,original))<=B
            self.assertEqual(feasible,sum(ci*bi for ci,bi in zip(m['costs'],reduced))<=m['budget'])
        q=core_qubo(m);z=1-2*q['bits']
        np.testing.assert_allclose(z@q['h']+np.einsum('bi,ij,bj->b',z,q['J'],z)+q['C'],q['E'],atol=1e-9)
        self.assertAlmostEqual(q['E'].min(),-float(exact(c,v,B,fixed)[1]))

    def test_presolve_total_zero_qubits(self):
        c=[2,3];v=[F(2),F(1)];r=presolve(c,v,5)
        solved=solve_core_qaoa(c,v,5,r['fixed'],shots=32)
        self.assertEqual(solved['qubits'],0);self.assertFalse(solved['execucao_quantica'])
        self.assertEqual(solved['bits'],[1,1])

    def test_pipeline_nucleo_com_ruido(self):
        # Empate: presolve não fixa nenhuma decisão ótima opcional.
        c=[1,1,1];v=[F(1,2)]*3;r=presolve(c,v,1)
        self.assertEqual(r['fixed'],{})
        sol=solve_core_qaoa(c,v,1,r['fixed'],shots=64,maxiter=20,restarts=1,noise=True)
        self.assertAlmostEqual(sum(x['prob'] for x in sol['distribuicao']),1.)
        self.assertEqual(sum(x['shots'] for x in sol['distribuicao_ruido']),64)

    def test_substituicoes_ising_equivalem_para_todos_estados(self):
        rng=np.random.default_rng(41);h=rng.normal(size=5);J=np.triu(rng.normal(size=(5,5)),1);C=.32
        for index in range(5):
            for other in [None,*[j for j in range(5) if j!=index]]:
                for sign in (-1,1):
                    nh,nJ,nC=eliminate(h,J,C,index,sign,other);keep=[i for i in range(5) if i!=index]
                    for reduced in product((-1,1),repeat=4):
                        z=np.zeros(5);z[keep]=reduced;z[index]=sign if other is None else sign*z[other]
                        original=z@h+z@J@z+C;reduced=np.array(reduced)
                        self.assertAlmostEqual(original,reduced@nh+reduced@nJ@reduced+nC)

    def test_recursivo_reconstroi_portfolio_e_qiro_preserva_orcamento(self):
        c=[1,2,3,2];v=[F(3,10),F(6,10),F(7,10),F(4,10)];B=4
        r=rqaoa(c,v,B,shots=64,maxiter=20,restarts=1,cutoff=3)
        self.assertEqual(len(r['bits']),4);self.assertFalse(r['preserva_otimo'])
        q=core_qubo(core(c,v,B));z=np.array(r['estado_Ising_reconstruido'])
        self.assertAlmostEqual(z@q['h']+z@q['J']@z+q['C'],r['energia_residual'])
        r=qiro(c,v,B,shots=64,maxiter=20,restarts=1,cutoff=2,backtracking=1)
        self.assertLessEqual(sum(ci*bi for ci,bi in zip(c,r['bits'])),B)
        self.assertFalse(r['preserva_otimo'])

    def test_agrupamento_marcado_como_heuristico(self):
        d=json.loads(Path(__file__).with_name('instancia.json').read_text())
        r=cluster_representatives(d,3)
        self.assertFalse(r['preserva_otimo']);self.assertEqual(len(r['representantes']),3)
        flattened=[i for g in r['grupos'] for i in g];self.assertEqual(sorted(flattened),list(range(7)))

if __name__=='__main__':unittest.main()
