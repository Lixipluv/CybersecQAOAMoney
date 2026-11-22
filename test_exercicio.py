import argparse
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from exercicio import (familia, instancia, codificar, interp, circuito_didatico, metricas, beneficio,
                       executar, calibrar_desbalanceada, LIMIAR)
from projeto import circuito, estado

CAMPOS_26 = ['ods', 'problem', 'instance', 'number_of_variables', 'number_of_qubits', 'penalty_lambda', 'p',
             'gamma', 'beta', 'optimizer', 'optimizer_iterations', 'shots', 'backend', 'seed', 'optimization_level']


class TestExercicio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = json.loads(Path(__file__).with_name('instancia.json').read_text())
        cls.ctrl = familia(cls.base, 10)

    def test_familia_aninhada_e_deterministica(self):
        self.assertEqual(self.ctrl[:7], self.base['controles'])
        self.assertEqual(familia(self.base, 10), self.ctrl)
        self.assertEqual(familia(self.base, 12)[:10], self.ctrl)
        for c in self.ctrl[7:]:
            self.assertTrue(c['sintetico']); self.assertIn(c['custo'], (2, 3, 4))
            self.assertGreaterEqual(sum(c['taticas']), 2); self.assertGreaterEqual(sum(c['ativos']), 1)

    def test_folga_qubo_ising_equivalentes(self):
        for k in (3, 7, 9):
            esp = codificar(instancia(self.base, self.ctrl, k), 'folga')
            self.assertLess(esp['erro_ising'], 1e-8); self.assertLess(esp['erro_qubo'], 1e-8)
            self.assertTrue(esp['ground_e_otimo'])
            self.assertAlmostEqual(-esp['E'].min(), esp['Vstar'])
            self.assertEqual(esp['n'], k + esp['B'].bit_length())

    def test_desbalanceada_sem_folga(self):
        esp = codificar(instancia(self.base, self.ctrl, 7), 'desbalanceada', (0.3, 0.1))
        self.assertEqual(esp['n'], 7)
        self.assertLess(esp['erro_ising'], 1e-8); self.assertLess(esp['erro_qubo'], 1e-8)
        h = esp['B'] - esp['custo']
        np.testing.assert_allclose(esp['E'], -esp['V'] - 0.3*h + 0.1*h*h)

    def test_calibracao_usa_so_instancias_de_calibracao(self):
        lam, info = calibrar_desbalanceada(self.base, self.ctrl, [3, 4])
        self.assertEqual(info['instancias_calibracao'], [3, 4])
        for k in (3, 4):
            self.assertTrue(codificar(instancia(self.base, self.ctrl, k), 'desbalanceada', lam)['ground_e_otimo'])

    def test_interp(self):
        np.testing.assert_allclose(interp(np.array([0.7])), [0.7, 0.7])
        np.testing.assert_allclose(interp(np.array([1., 2.])), [1., 1.5, 2.])

    def test_circuito_didatico_equivale(self):
        from qiskit.quantum_info import Operator
        esp = codificar(instancia(self.base, self.ctrl, 3), 'folga')
        params = np.array([0.4, 1.1, 0.3, 0.9])
        escala = max(float(np.max(np.abs(esp['E']))), 1.)
        a = circuito(esp['h'], esp['J'], esp['C'], params, esp['E'])
        b = circuito_didatico(esp['h'], esp['J'], params, escala)
        b.remove_final_measurements(inplace=True)
        b.data = [i for i in b.data if i.operation.name != 'barrier']
        self.assertTrue(Operator(a).equiv(Operator(b)))

    def test_metricas_limites(self):
        esp = codificar(instancia(self.base, self.ctrl, 4), 'folga')
        unif = np.full(2**esp['n'], 1/2**esp['n'])
        m = metricas(esp, prob=unif, counts=np.random.default_rng(1).multinomial(4096, unif))
        self.assertAlmostEqual(m['ganho95_vs_uniforme'], 1.)
        self.assertFalse(m['util'])
        pico = esp['opt'] & esp['ground']; prob = pico/pico.sum()
        m = metricas(esp, prob=prob, counts=np.random.default_rng(1).multinomial(100, prob))
        self.assertEqual(m['Popt'], 1.); self.assertEqual(m['R_app'], 1.)
        self.assertTrue(m['util']); self.assertTrue(m['menor_energia_e_mais_frequente'])

    def test_beneficio_referencia(self):
        esp = codificar(instancia(self.base, self.ctrl, 7), 'folga')
        cen = beneficio(esp, [])
        self.assertEqual(cen['referencia_ordem_CIS']['delta_valor_pct_vs_ref'], 0.)
        self.assertAlmostEqual(cen['otimo_classico']['valor'], esp['Vstar'])
        self.assertGreaterEqual(cen['otimo_classico']['valor'], cen['gulosa_valor_custo']['valor'])
        self.assertLessEqual(cen['referencia_ordem_CIS']['custo'], esp['B'])

    def test_execucao_com_fake_kingston_registra_reprodutibilidade(self):
        esp = codificar(instancia(self.base, self.ctrl, 3), 'folga')
        args = argparse.Namespace(maxiter=10, restarts=1, shots=128, nivel=1)
        with tempfile.TemporaryDirectory() as tmp:
            r, params, tc, med = executar(esp, 1, 42, args, None, Path(tmp), e2=True)
            salvo = json.loads((Path(tmp)/'resultado.json').read_text())
            for env in ('E1', 'E2'):
                for campo in CAMPOS_26:
                    self.assertIn(campo, salvo[env]['reprodutibilidade'])
            self.assertEqual(salvo['E2']['reprodutibilidade']['optimization_level'], 1)
            self.assertEqual(salvo['E2']['shots']['shots'], 128)
            self.assertTrue((Path(tmp)/'medicoes_E2.csv').exists())
            self.assertTrue((Path(tmp)/'circuito_fake_kingston.qasm').exists())
            self.assertGreater(r['recursos']['fake_kingston']['gates_2q'], 0)
            self.assertTrue(set(r['recursos']['fake_kingston']['por_tipo']) <= {'rz', 'sx', 'x', 'cz'})


if __name__ == '__main__':
    unittest.main()
