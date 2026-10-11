import json
from pathlib import Path
import unittest
import numpy as np
from mixer_factivel import prepare,state,build

class TestMixer(unittest.TestCase):
    def test_estado_e_circuito_preservam_factivel(self):
        from qiskit.quantum_info import Statevector
        d=json.loads(Path(__file__).with_name('instancia.json').read_text())
        d['controles']=d['controles'][:4];d['orcamento']=4
        R,b,cost,V,feasible,initial=prepare(d);scale=max(float(R.sum()),1.)
        for p in (1,2,3):
            params=np.random.default_rng(p).uniform(0,7,2*p)
            psi=state(params,-V,initial,scale)
            self.assertAlmostEqual(float(np.sum(abs(psi)**2)),1.)
            self.assertLess(float(np.sum(abs(psi[~feasible])**2)),1e-20)
            q=Statevector.from_instruction(build(R,initial,params,scale)).data
            np.testing.assert_allclose(psi,q,atol=1e-9)

if __name__=='__main__':unittest.main()
