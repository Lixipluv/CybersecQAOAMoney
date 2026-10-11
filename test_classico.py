import copy
import json
from pathlib import Path
import unittest
import numpy as np
from classico import resolver
from projeto import tabela

class TestClassico(unittest.TestCase):
    def test_dp_equivale_exaustivo(self):
        source=json.loads(Path(__file__).with_name('instancia.json').read_text())
        rng=np.random.default_rng(93)
        for _ in range(20):
            d=copy.deepcopy(source)
            for c in d['controles']:c['custo']=int(rng.integers(1,6))
            d['orcamento']=int(rng.integers(0,15))
            *_,bits,cost,V,E=tabela(d)
            result=resolver(d)
            self.assertAlmostEqual(result['valor'],V[cost<=d['orcamento']].max())
            self.assertLessEqual(result['custo'],d['orcamento'])

if __name__=='__main__':unittest.main()
