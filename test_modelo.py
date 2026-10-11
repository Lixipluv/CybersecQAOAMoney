import json
from pathlib import Path
import unittest
import tempfile
import numpy as np
from projeto import tabela, estado, circuito, executar

class TestModelo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads(Path(__file__).with_name('instancia.json').read_text())

    def test_energias_e_solucao_original(self):
        R,a,rho,Q,h,J,C,b,costs,V,E=tabela(self.d)
        z=1-2*b
        np.testing.assert_allclose(np.einsum('bi,ij,bj->b',b,Q,b)+rho*self.d['orcamento']**2,E,atol=1e-9)
        np.testing.assert_allclose(z@h+np.einsum('bi,ij,bj->b',z,J,z)+C,E,atol=1e-9)
        ground=np.isclose(E,E.min())
        self.assertTrue(np.all(costs[ground]<=self.d['orcamento']))
        self.assertAlmostEqual(-E.min(),V[costs<=self.d['orcamento']].max())

    def test_statevector_equivale_qiskit(self):
        from qiskit.quantum_info import Statevector
        R,a,rho,Q,h,J,C,b,costs,V,E=tabela(self.d)
        for p in [1,2,3]:
            params=np.random.default_rng(p).uniform(0,2,2*p)
            ours=estado(params,E,len(a))
            qiskit=Statevector.from_instruction(circuito(h,J,C,params,E)).data
            np.testing.assert_allclose(ours,qiskit,atol=1e-10)
            self.assertAlmostEqual(float(np.vdot(ours,ours).real),1.)

    def test_exportacao_completa(self):
        d=dict(self.d); d['controles']=d['controles'][:3]; d['orcamento']=3
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            r=executar(d,1,42,64,20,1,out)
            saved=json.loads((out/'resultado.json').read_text())
            self.assertEqual(saved['shots'],64)
            self.assertTrue((out/'medicoes.csv').exists())
            self.assertTrue(0 <= r['Pfeas'] <= 1)

    def test_folga_orcamento_zero(self):
        d=dict(self.d);d['orcamento']=0
        *_,E=tabela(d)
        self.assertAlmostEqual(E.min(),0.)

    def test_leitura_transpilada_preserva_bits(self):
        from qiskit import QuantumCircuit, transpile
        from qiskit.transpiler import CouplingMap
        from qiskit_aer import AerSimulator
        qc=QuantumCircuit(5);qc.x(0);qc.x(3);qc.cx(0,4);qc.measure_all()
        tc=transpile(qc,basis_gates=['rz','sx','x','cx'],coupling_map=CouplingMap.from_line(5),seed_transpiler=42)
        counts=AerSimulator().run(tc,shots=32).result().get_counts()
        self.assertEqual(counts,{'11001':32})

if __name__=='__main__':unittest.main()
