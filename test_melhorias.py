import copy
import json
import tempfile
from pathlib import Path
import unittest
import numpy as np
from projeto import tabela, estado, circuito
from melhorias import cvar, greedy, run

class TestMelhorias(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads(Path(__file__).with_name('instancia.json').read_text())

    def test_cvar_massa_fracionaria(self):
        prob=np.array([.1,.4,.5]); E=np.array([-2.,1.,10.])
        self.assertAlmostEqual(cvar(prob,E,.25),(-.2+.15)/.25)
        self.assertAlmostEqual(cvar(prob,E,1),prob@E)
        self.assertAlmostEqual(cvar(prob,E,.05),-2.)

    def test_penalidade_menor_preserva_otimo(self):
        from projeto import valores
        rng=np.random.default_rng(17)
        for _ in range(30):
            d=copy.deepcopy(self.d);d['controles']=d['controles'][:4]
            for c in d['controles']:c['custo']=int(rng.integers(1,5))
            d['orcamento']=int(rng.integers(0,12))
            rho=float(valores(d).max()+.05)
            R,a,rho,Q,h,J,C,b,costs,V,E=tabela(d,rho)
            ground=np.isclose(E,E.min(),atol=1e-9,rtol=0)
            self.assertTrue(np.all(costs[ground]<=d['orcamento']))
            self.assertAlmostEqual(-E.min(),V[costs<=d['orcamento']].max())
            z=1-2*b
            np.testing.assert_allclose(z@h+np.einsum('bi,ij,bj->b',z,J,z)+C,E,atol=1e-9)

    def test_warm_circuito_e_estado(self):
        from qiskit.quantum_info import Statevector
        R,a,rho,Q,h,J,C,b,c,V,E=tabela(self.d)
        n=len(a);rng=np.random.default_rng(7)
        for p in (1,2,3):
            warm=rng.uniform(.1,.9,n);params=rng.uniform(0,3,2*p)
            ours=estado(params,E,n,scale=17.,warm=warm)
            qs=Statevector.from_instruction(circuito(h,J,C,params,E,scale=17.,warm=warm)).data
            np.testing.assert_allclose(ours,qs,atol=1e-10)
        # Mixer alinhado deixa a preparação inicial invariante até fase global.
        base=estado(np.zeros(2),E,n,17.,warm)
        mixed=estado(np.array([0.,.7]),E,n,17.,warm)
        self.assertAlmostEqual(float(abs(np.vdot(base,mixed))),1.)

    def test_greedy_factivel_e_folga(self):
        from projeto import valores
        R=valores(self.d);g=greedy(self.d,R)
        k=len(R);cost=sum(c['custo']*g[i] for i,c in enumerate(self.d['controles']))
        slack=sum(2**m*g[k+m] for m in range(len(g)-k))
        self.assertLessEqual(cost,self.d['orcamento'])
        self.assertEqual(cost+slack,self.d['orcamento'])

    def test_pipeline_ruido(self):
        d=copy.deepcopy(self.d);d['controles']=d['controles'][:3];d['orcamento']=3
        with tempfile.TemporaryDirectory() as tmp:
            r=run(d,'warm_cvar',1,42,64,20,1,Path(tmp)/'run',noise=True)
            self.assertTrue(0<=r['Popt']<=1)
            self.assertTrue(0<=r['ruido']['Pfeas']<=1)
            self.assertTrue((Path(tmp)/'run'/'resultado.json').exists())

if __name__=='__main__':unittest.main()
