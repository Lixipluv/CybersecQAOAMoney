"""Protocolo completo do exercício OQA (seções 5–26 do enunciado).

E1: simulação ideal (statevector NumPy validado contra Qiskit; probabilidades exatas e shots).
E2: FakeKingston (IBM Heron r2): transpilação ISA e ruído de calibração (portas, leitura,
    relaxação térmica) via AerSimulator.from_backend.

Instâncias: prefixos da base CIS (3–7 controles, dados de instancia.json) e extensão
sintética aninhada (controles S8, S9, ...) com atributos sorteados nas mesmas faixas da base.

Codificações:
- folga: QUBO padrão do exercício, restrição de orçamento via folga binária e penalidade rho;
- desbalanceada: penalização desbalanceada sem folga, E=-V-l1*h+l2*h^2, h=B-sum(c x).
  (l1, l2) são calibrados apenas nas instâncias de calibração e depois fixados.

O ótimo clássico nunca participa da otimização variacional; só da avaliação.
"""
import argparse
import copy
import csv
import json
import math
import platform
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from projeto import tabela, estado, circuito

ROOT = Path(__file__).resolve().parent
ODS = {'numero': 16, 'nome': 'Paz, Justiça e Instituições Eficazes',
       'metas': ['16.6: instituições eficazes, responsáveis e transparentes',
                 '16.10: acesso público à informação e proteção das liberdades fundamentais']}
PROBLEMA = 'Alocação de recursos de cibersegurança: seleção de controles CIS sob orçamento (mochila 0-1)'
LIMIAR = 0.95            # qualidade mínima: 95% do valor ótimo
OTIMIZADOR = 'COBYLA'
BACKEND_E1 = 'statevector_numpy_ideal'
BACKEND_E2 = 'aer_from_backend_fake_kingston'


# ----------------------------------------------------------------------------- instâncias
def familia(base, maximo, seed=2026):
    """Base real + controles sintéticos determinísticos; prefixos ficam aninhados."""
    rng = np.random.default_rng(seed)
    controles = copy.deepcopy(base['controles'])
    nt, na = len(base['relevancia_taticas']), len(base['ativos'])
    for i in range(len(controles), maximo):
        c = copy.deepcopy(base['controles'][i % len(base['controles'])])
        c['nome'] = f'S{i+1}: controle sintético'
        c['custo'] = int(rng.integers(2, 5))
        c['eficacia'] = round(float(rng.uniform(.65, .9)), 3)
        c['adesao'] = round(float(rng.uniform(.7, .9)), 3)
        c['prioridade'] = round(float(rng.uniform(.8, .95)), 3)
        c['impacto'] = round(float(rng.uniform(.15, .4)), 3)
        tat = np.zeros(nt, int); tat[rng.choice(nt, int(rng.integers(2, 4)), replace=False)] = 1
        atv = np.zeros(na, int); atv[rng.choice(na, int(rng.integers(1, na+1)), replace=False)] = 1
        c['taticas'], c['ativos'], c['sintetico'] = tat.tolist(), atv.tolist(), True
        controles.append(c)
    return controles


def instancia(base, controles, k):
    d = copy.deepcopy(base)
    d['controles'] = copy.deepcopy(controles[:k])
    d['orcamento'] = max(1, int(sum(c['custo'] for c in d['controles']) * base['fracao_orcamento']))
    return d


def nome_instancia(k):
    return f'cis_prefixo_{k}' if k <= 7 else f'cis7_mais_sinteticos_{k}'


# ----------------------------------------------------------------------------- codificações
def qubo_para_ising(Q, c0):
    """Q triangular superior (diagonal = termos lineares); x=(1-z)/2."""
    lin = np.diag(Q).copy()
    J = np.triu(Q, 1)/4
    h = -lin/2 - (J.sum(axis=0)+J.sum(axis=1))
    C = c0 + lin.sum()/2 + J.sum()
    return h, J, C


def _bits(n):
    idx = np.arange(2**n)
    return ((idx[:, None] >> np.arange(n)) & 1).astype(float)


def codificar(d, tipo='folga', lambdas=None):
    """Monta QUBO, Ising e o espaço completo de estados com todas as máscaras de avaliação."""
    if tipo == 'folga':
        R, a, rho, Q, h, J, C, bits, costs, V, E = tabela(d)
        c0 = rho*d['orcamento']**2
        penal = {'rho': rho}
    elif tipo == 'desbalanceada':
        from projeto import valores
        R = valores(d); a = np.array([c['custo'] for c in d['controles']], float)
        B = d['orcamento']; l1, l2 = lambdas
        # E(x) = -R.x - l1 (B - a.x) + l2 (B - a.x)^2
        lin = -R + l1*a + l2*(a*a - 2*B*a)
        Q = np.diag(lin)
        for i in range(len(a)):
            for j in range(i+1, len(a)):
                Q[i, j] = 2*l2*a[i]*a[j]
        c0 = -l1*B + l2*B*B
        h, J, C = qubo_para_ising(Q, c0)
        bits = _bits(len(a)); costs = bits @ a; V = bits @ R
        E = -V - l1*(B-costs) + l2*(B-costs)**2
        penal = {'lambda1': l1, 'lambda2': l2}
    else:
        raise ValueError(tipo)
    k, n, B = len(R), len(h), d['orcamento']
    feas = costs <= B
    Vs = float(V[feas].max())
    opt = feas & np.isclose(V, Vs, atol=1e-10, rtol=0)
    good = feas & (V >= LIMIAR*Vs - 1e-12)
    ground = np.isclose(E, E.min(), atol=1e-9, rtol=0)
    q = np.array(d['relevancia_taticas'], float)
    T = np.array([c['taticas'] for c in d['controles']], float)
    cov = (np.minimum(1, bits[:, :k] @ T) @ q)/q.sum()
    acc = feas & (V >= d['valor_minimo']) & (cov >= d['cobertura_minima'])
    z = 1-2*bits
    e_ising = (z @ h) + np.einsum('bi,bi->b', z @ J, z) + C
    e_qubo = np.einsum('bi,bi->b', bits @ Q, bits) + c0
    return {'tipo': tipo, 'd': d, 'k': k, 'n': n, 'B': B, 'R': R, 'custos_ctrl': np.array([c['custo'] for c in d['controles']], float),
            'penalidade': penal, 'Q': Q, 'c0': c0, 'h': h, 'J': J, 'C': C, 'bits': bits, 'custo': costs, 'V': V, 'E': E,
            'feas': feas, 'opt': opt, 'good': good, 'ground': ground, 'cov': cov, 'acc': acc, 'Vstar': Vs,
            'erro_ising': float(np.max(np.abs(e_ising-E))), 'erro_qubo': float(np.max(np.abs(e_qubo-E))),
            'ground_e_otimo': bool(np.all(opt[ground])),
            'uniforme': {'Popt': float(opt.mean()), 'Pfeas': float(feas.mean()), 'P95': float(good.mean()),
                         'Pqubo': float(ground.mean()), 'Paceitavel': float(acc.mean())}}


def caracteristicas(esp):
    """Seção 12.1: tamanho e estrutura do modelo."""
    Q = esp['Q']; n = esp['n']
    quad = int(np.count_nonzero(np.triu(Q, 1)))
    return {'controles': esp['k'], 'variaveis_decisao': esp['k'], 'variaveis_auxiliares': n-esp['k'],
            'variaveis_total': n, 'restricoes': 1, 'orcamento': esp['B'],
            'termos_lineares': int(np.count_nonzero(np.abs(esp['h']) > 1e-12)),
            'termos_quadraticos': quad, 'interacoes_ZZ': int(np.count_nonzero(np.abs(esp['J']) > 1e-12)),
            'densidade_qubo': quad/(n*(n-1)/2) if n > 1 else 0., 'estados': 2**n,
            'portfolios_factiveis': int(len({tuple(b) for b in esp['bits'][esp['feas'], :esp['k']].astype(int)})),
            'portfolios_totais': 2**esp['k'], 'escala_max_abs_E': float(np.max(np.abs(esp['E']))),
            'razao_escala_valor_otimo': float(np.max(np.abs(esp['E']))/esp['Vstar']),
            'memoria_statevector_MB': 2**n*16/2**20, **esp['penalidade']}


def calibrar_desbalanceada(base, controles, ks):
    """Escolhe (l1,l2) só nas instâncias de calibração: estado fundamental = ótimo em todas e
    maior gap relativo mínimo (E_1º não ótimo - E_ótimo)/(max E - min E). Sem usar instâncias de teste."""
    grade = [(round(float(l1), 3), round(float(l2), 3)) for l1 in np.linspace(.05, 1., 20) for l2 in np.linspace(.01, .5, 50)]
    melhor, registro = None, []
    for lam in grade:
        gaps = []
        for k in ks:
            esp = codificar(instancia(base, controles, k), 'desbalanceada', lam)
            if not esp['ground_e_otimo']:
                gaps = None; break
            E = esp['E']; nao = E[~esp['opt']]
            gaps.append((nao.min()-E.min())/(E.max()-E.min()))
        if gaps is not None:
            registro.append({'lambda1': lam[0], 'lambda2': lam[1], 'gap_min': min(gaps)})
            if melhor is None or min(gaps) > melhor[1]:
                melhor = (lam, min(gaps))
    if melhor is None:
        raise RuntimeError('Nenhum (l1,l2) da grade preserva o ótimo na calibração.')
    return melhor[0], {'instancias_calibracao': ks, 'pontos_validos': len(registro), 'pontos_grade': len(grade),
                       'escolhido': {'lambda1': melhor[0][0], 'lambda2': melhor[0][1], 'gap_relativo_min': melhor[1]}}


# ----------------------------------------------------------------------------- QAOA
def interp(v):
    """Inicialização INTERP (Zhou et al., 2020): parâmetros de p para p+1."""
    p = len(v); w = np.r_[0, v, 0]
    return np.array([(i/p)*w[i] + ((p-i)/p)*w[i+1] for i in range(p+1)])


def otimizar(E, n, p, seed, maxiter, restarts, init=None):
    rng = np.random.default_rng(seed)
    escala = max(float(np.max(np.abs(E))), 1.)
    inicios = [np.asarray(init, float)] if init is not None else []
    while len(inicios) < restarts:
        inicios.append(np.r_[rng.uniform(0, 2*np.pi, p), rng.uniform(0, np.pi, p)])
    runs = []
    for i, x0 in enumerate(inicios):
        trace = []
        def f(x):
            e = float(np.abs(estado(x, E, n, escala))**2 @ E)
            trace.append(e); return e/escala
        t = time.perf_counter()
        res = minimize(f, x0, method=OTIMIZADOR, options={'maxiter': maxiter, 'rhobeg': 1.0, 'tol': 1e-5})
        runs.append({'origem': 'interp' if (init is not None and i == 0) else 'aleatorio',
                     'iniciais': x0.tolist(), 'finais': res.x.tolist(), 'avaliacoes': int(res.nfev),
                     'energia_inicial': trace[0], 'energia_final': float(res.fun)*escala,
                     'convergiu': bool(res.success), 'mensagem': str(res.message),
                     'tempo_s': time.perf_counter()-t, 'traco': trace})
    best = min(range(len(runs)), key=lambda i: runs[i]['energia_final'])
    return np.array(runs[best]['finais']), runs, best, escala


def circuito_didatico(h, J, params, escala):
    """Mesmo circuito de projeto.circuito, com barreiras rotuladas para a figura."""
    from qiskit import QuantumCircuit
    n, p = len(h), len(params)//2
    qc = QuantumCircuit(n)
    qc.h(range(n)); qc.barrier(label='|+>^n')
    for l, (g, b) in enumerate(zip(params[:p], params[p:]), 1):
        for i in range(n):
            if abs(h[i]) > 1e-12: qc.rz(2*g*h[i]/escala, i)
        for i in range(n):
            for j in range(i+1, n):
                if abs(J[i, j]) > 1e-12: qc.rzz(2*g*J[i, j]/escala, i, j)
        qc.barrier(label=f'U_C(γ{l})')
        qc.rx(2*b, range(n)); qc.barrier(label=f'U_M(β{l})')
    qc.measure_all(add_bits=True)
    return qc


_BACKEND = {}
def fake_kingston():
    if 'b' not in _BACKEND:
        from qiskit_ibm_runtime.fake_provider import FakeKingston
        _BACKEND['b'] = FakeKingston()
    return _BACKEND['b']


def transpilar(qc, nivel, seed, roteado=True):
    from qiskit import transpile
    if roteado:
        from qiskit.transpiler import generate_preset_pass_manager
        return generate_preset_pass_manager(optimization_level=nivel, backend=fake_kingston(),
                                            seed_transpiler=seed).run(qc)
    # Mesma base nativa (cz, rz, sx, x), conectividade total: isola o custo do roteamento.
    return transpile(qc, basis_gates=['cz', 'rz', 'sx', 'x'], optimization_level=nivel, seed_transpiler=seed)


def recursos(qc):
    ops = [i for i in qc.data if i.operation.name not in ('measure', 'barrier', 'delay')]
    ativos = {qc.find_bit(q).index for i in qc.data if i.operation.name != 'barrier' for q in i.qubits}
    tipos = {}
    for i in ops: tipos[i.operation.name] = tipos.get(i.operation.name, 0)+1
    return {'qubits': len(ativos), 'profundidade': qc.depth(lambda i: i.operation.name not in ('measure', 'barrier')),
            'gates': len(ops), 'gates_2q': sum(len(i.qubits) == 2 for i in ops), 'por_tipo': tipos}


def simular_ruido(tc, shots, seed, nm=None):
    from qiskit_aer import AerSimulator
    if nm is None:
        if 'sim' not in _BACKEND: _BACKEND['sim'] = AerSimulator.from_backend(fake_kingston())
        sim = _BACKEND['sim']
    else:
        sim = AerSimulator(noise_model=nm)
    raw = sim.run(tc, shots=shots, seed_simulator=seed).result().get_counts()
    counts = None
    for b, c in raw.items():
        b = b.replace(' ', '')
        if counts is None: counts = np.zeros(2**len(b), int)
        counts[int(b, 2)] += c
    return counts


def wilson_inf(x, n, z=1.96):
    if n == 0: return 0.
    ph = x/n; den = 1+z*z/n
    return max(0., (ph + z*z/(2*n) - z*math.sqrt(ph*(1-ph)/n + z*z/(4*n*n)))/den)


def decodificar(esp, i):
    if i is None: return None
    k, n = esp['k'], esp['n']
    return {'bitstring': format(int(i), f'0{n}b'),
            'decisao': ''.join(str(int(b)) for b in esp['bits'][i, :k]),
            'controles': [esp['d']['controles'][j]['nome'] for j in range(k) if esp['bits'][i, j]],
            'custo': float(esp['custo'][i]), 'valor': float(esp['V'][i]), 'cobertura': float(esp['cov'][i]),
            'factivel': bool(esp['feas'][i]), 'otimo': bool(esp['opt'][i]), 'aceitavel': bool(esp['acc'][i]),
            'energia': float(esp['E'][i])}


def metricas(esp, prob=None, counts=None):
    """Probabilidades exatas (prob) ou frequências (counts); amostras sempre via counts."""
    shots = int(counts.sum())
    P = prob if prob is not None else counts/shots
    m = {k: float(P[esp[mk]].sum()) for k, mk in
         [('Popt', 'opt'), ('Pfeas', 'feas'), ('Pqubo', 'ground'), ('P95', 'good'), ('Paceitavel', 'acc')]}
    m['energia_esperada'] = float(P @ esp['E'])
    amostr = np.flatnonzero(counts)
    melhor = int(amostr[np.argmin(esp['E'][amostr])])
    fa = amostr[esp['feas'][amostr]]
    melhor_f = int(fa[np.argmax(esp['V'][fa])]) if len(fa) else None
    moda = int(np.argmax(counts))
    m['R_app'] = float(esp['V'][melhor_f]/esp['Vstar']) if melhor_f is not None else 0.
    m['otimo_observado'] = bool(esp['opt'][amostr].any())
    m['E_QAOA'] = float(esp['E'][melhor])
    m['melhor_energia'] = decodificar(esp, melhor)
    m['melhor_factivel'] = decodificar(esp, melhor_f)
    m['mais_frequente'] = decodificar(esp, moda)
    m['menor_energia_e_mais_frequente'] = bool(np.isclose(esp['E'][melhor], esp['E'][moda]))
    u = esp['uniforme']['P95']
    m['ganho95_vs_uniforme'] = m['P95']/u if u > 0 else float('nan')
    m['ganho_opt_vs_uniforme'] = m['Popt']/esp['uniforme']['Popt']
    lb = m['P95'] if prob is not None else wilson_inf(int(counts[esp['good']].sum()), shots)
    m['P95_limite_inferior_95'] = lb
    m['shots_para_95pct_de_obter_95'] = (math.ceil(math.log(.05)/math.log1p(-m['P95'])) if 0 < m['P95'] < 1
                                         else (1 if m['P95'] >= 1 else None))
    m['decisao_util'] = bool(m['R_app'] >= LIMIAR)                 # nível A: decisão útil obtida
    m['util'] = bool(m['decisao_util'] and lb > u)                  # nível B: útil e atribuível ao QAOA
    m['fonte'] = 'probabilidade_exata' if prob is not None else 'frequencia_shots'
    m['shots'] = shots
    return m


def reprodutibilidade(esp, p, params, runs, shots, backend, seed, nivel):
    p_ = len(params)//2
    return {'ods': f"ODS {ODS['numero']} – {ODS['nome']}", 'problem': PROBLEMA, 'instance': nome_instancia(esp['k']),
            'encoding': esp['tipo'], 'number_of_variables': esp['n'], 'number_of_decision_variables': esp['k'],
            'number_of_auxiliary_variables': esp['n']-esp['k'], 'number_of_qubits': esp['n'],
            'penalty_lambda': esp['penalidade'], 'p': p, 'gamma': list(map(float, params[:p_])),
            'beta': list(map(float, params[p_:])), 'optimizer': OTIMIZADOR,
            'optimizer_iterations': int(sum(r['avaliacoes'] for r in runs)), 'shots': shots,
            'backend': backend, 'seed': seed, 'optimization_level': nivel}


def executar(esp, p, seed, args, init, out, e2=True):
    t0 = time.perf_counter()
    n, E = esp['n'], esp['E']
    params, runs, best, escala = otimizar(E, n, p, seed, args.maxiter*p, args.restarts, init)
    t_opt = sum(r['tempo_s'] for r in runs)
    t = time.perf_counter()
    prob = np.abs(estado(params, E, n, escala))**2
    counts1 = np.random.default_rng(seed+10_000).multinomial(args.shots, prob)
    t_e1 = time.perf_counter()-t
    qc = circuito(esp['h'], esp['J'], esp['C'], params, E)
    med = qc.copy(); med.measure_all()
    t = time.perf_counter(); tc = transpilar(med, args.nivel, seed); t_tr = time.perf_counter()-t
    tb = transpilar(med, args.nivel, seed, roteado=False)
    r = {'ods': ODS, 'problema': PROBLEMA, 'instancia': nome_instancia(esp['k']), 'codificacao': esp['tipo'],
         'k': esp['k'], 'n': esp['n'], 'p': p, 'seed': seed,
         'E_otima': float(E.min()), 'valor_otimo': esp['Vstar'], 'escala_hamiltoniano': escala,
         'classico': decodificar(esp, int(np.flatnonzero(esp['opt'] & esp['ground'])[0]) if (esp['opt'] & esp['ground']).any()
                                 else int(np.flatnonzero(esp['opt'])[0])),
         'otimizacao': {'escolhido': best, 'inicios': [{k: v for k, v in x.items() if k != 'traco'} for x in runs],
                        'tempo_s': t_opt, 'avaliacoes_total': int(sum(x['avaliacoes'] for x in runs))},
         'recursos': {'logico': recursos(med), 'base_nativa_sem_roteamento': recursos(tb), 'fake_kingston': recursos(tc)},
         'tempo_transpilacao_s': t_tr,
         'E1': {'reprodutibilidade': reprodutibilidade(esp, p, params, runs, args.shots, BACKEND_E1, seed, None),
                'exata': metricas(esp, prob=prob, counts=counts1), 'shots': metricas(esp, counts=counts1),
                'tempo_execucao_s': t_e1}}
    counts2 = None
    if e2:
        t = time.perf_counter(); counts2 = simular_ruido(tc, args.shots, seed); t_e2 = time.perf_counter()-t
        r['E2'] = {'reprodutibilidade': reprodutibilidade(esp, p, params, runs, args.shots, BACKEND_E2, seed, args.nivel),
                   'shots': metricas(esp, counts=counts2), 'tempo_execucao_s': t_e2}
        r['delta_ruido'] = {k: r['E1']['exata'][k]-r['E2']['shots'][k] for k in ('Popt', 'Pfeas', 'P95', 'Pqubo')}
    r['tempo_total_s'] = time.perf_counter()-t0
    out.mkdir(parents=True, exist_ok=True)
    (out/'resultado.json').write_text(json.dumps(r, indent=2, ensure_ascii=False), encoding='utf-8')
    from qiskit import qasm3
    (out/'circuito_fake_kingston.qasm').write_text(qasm3.dumps(tc), encoding='utf-8')
    for nome, cts in (('E1', counts1), ('E2', counts2)):
        if cts is None: continue
        with (out/f'medicoes_{nome}.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['bitstring', 'decisao', 'shots', 'energia', 'factivel', 'otimo', 'valor', 'custo', 'cobertura'])
            for i in np.flatnonzero(cts):
                dd = decodificar(esp, i)
                w.writerow([dd['bitstring'], dd['decisao'], int(cts[i]), dd['energia'], dd['factivel'], dd['otimo'], dd['valor'], dd['custo'], dd['cobertura']])
    figura_energia(runs, E.min(), esp, p, out/'energia_iteracao.png')
    figura_histograma(esp, counts1, counts2, p, out/'histograma.png')
    return r, params, tc, med


# ----------------------------------------------------------------------------- benefício ODS
def _carteira(esp, x):
    x = np.asarray(x, float); k = esp['k']
    T = np.array([c['taticas'] for c in esp['d']['controles']], float); q = np.array(esp['d']['relevancia_taticas'], float)
    return {'decisao': ''.join(str(int(v)) for v in x), 'custo': float(x @ esp['custos_ctrl']), 'valor': float(x @ esp['R']),
            'cobertura': float((np.minimum(1, x @ T) @ q)/q.sum()),
            'controles': [esp['d']['controles'][j]['nome'] for j in range(k) if x[j]]}


def primeiro_ajuste(esp, ordem):
    x = np.zeros(esp['k']); gasto = 0
    for j in ordem:
        if gasto + esp['custos_ctrl'][j] <= esp['B']:
            x[j] = 1; gasto += esp['custos_ctrl'][j]
    return x


def beneficio(esp, runs):
    """Seção 22: cenário de referência x cenário otimizado (valor de segurança e cobertura ATT&CK)."""
    k = esp['k']
    xs = _bits(k); custos = xs @ esp['custos_ctrl']; fe = custos <= esp['B']
    T = np.array([c['taticas'] for c in esp['d']['controles']], float); q = np.array(esp['d']['relevancia_taticas'], float)
    cen = {'referencia_ordem_CIS': _carteira(esp, primeiro_ajuste(esp, range(k))),
           'gulosa_valor_custo': _carteira(esp, primeiro_ajuste(esp, np.argsort(-esp['R']/esp['custos_ctrl'], kind='stable'))),
           'aleatoria_factivel_media': {'valor': float((xs[fe] @ esp['R']).mean()),
                                        'cobertura': float(((np.minimum(1, xs[fe] @ T) @ q)/q.sum()).mean()),
                                        'custo': float(custos[fe].mean())},
           'otimo_classico': _carteira(esp, esp['bits'][np.flatnonzero(esp['opt'])[0], :k])}
    for env in ('E1', 'E2'):
        rs = [r for r in runs if env in r]
        if not rs: continue
        mets = [r[env]['shots'] for r in rs]
        bf = [m['melhor_factivel'] for m in mets if m['melhor_factivel']]
        if bf:
            cen[f'QAOA_{env}_melhor_amostra_media'] = {'valor': float(np.mean([b['valor'] for b in bf])),
                                                     'cobertura': float(np.mean([b['cobertura'] for b in bf])),
                                                     'custo': float(np.mean([b['custo'] for b in bf]))}
        mf = [m['mais_frequente'] for m in mets]
        cen[f'QAOA_{env}_mais_frequente_media'] = {'valor': float(np.mean([b['valor'] if b['factivel'] else 0 for b in mf])),
                                                 'cobertura': float(np.mean([b['cobertura'] if b['factivel'] else 0 for b in mf])),
                                                 'custo': float(np.mean([b['custo'] for b in mf]))}
    ref = cen['referencia_ordem_CIS']
    for c in cen.values():
        c['delta_valor_pct_vs_ref'] = 100*(c['valor']-ref['valor'])/ref['valor']
        c['delta_cobertura_pp_vs_ref'] = 100*(c['cobertura']-ref['cobertura'])
        c['fracao_do_otimo'] = c['valor']/esp['Vstar']
    return cen


# ----------------------------------------------------------------------------- figuras
CORES = {'E1': '#2f6db5', 'E2': '#c4512b', 'unif': '#8a8a8a', 'desb': '#1c9275'}


def figura_energia(runs, emin, esp, p, path):
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for i, r in enumerate(runs):
        ax.plot(r['traco'], lw=1, label=f"início {i+1} ({r['origem']})")
    ax.axhline(emin, color='k', ls='--', lw=1, label='E* (ótimo QUBO)')
    ax.set(xlabel='Avaliação da função (COBYLA)', ylabel='⟨H_C⟩', title=f"{esp['tipo']}: {esp['k']} controles, {esp['n']} qubits, p={p}")
    ax.set_yscale('symlog', linthresh=max(1., abs(emin)))
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def figura_histograma(esp, c1, c2, p, path):
    """Linha 1: bitstrings medidos. Linha 2: agregado por portfólio (decisão), com referência uniforme."""
    envs = [('E1 ideal (shots)', c1)] + ([('E2 FakeKingston', c2)] if c2 is not None else [])
    k = esp['k']; dec = (esp['bits'][:, :k] @ (2**np.arange(k))).astype(int)
    xs = _bits(k); fe_d = (xs @ esp['custos_ctrl']) <= esp['B']
    opt_d = int(dec[np.flatnonzero(esp['opt'])[0]])
    fig, axes = plt.subplots(2, len(envs), figsize=(6.5*len(envs), 7), squeeze=False)
    for col, (nome, cts) in enumerate(envs):
        ax = axes[0, col]
        top = np.argsort(cts)[::-1][:20]; top = top[cts[top] > 0]
        cor = ['#1c9275' if esp['opt'][i] else '#2f6db5' if esp['feas'][i] else '#d0614c' for i in top]
        ax.bar(range(len(top)), cts[top], color=cor)
        ax.set_xticks(range(len(top)), [format(int(i), f"0{esp['n']}b") for i in top], rotation=90, fontsize=6)
        ax.set_title(f"{nome}: bitstrings (ótimo clássico observado {int(cts[esp['opt']].sum())}× de {int(cts.sum())})", fontsize=9)
        ax.set_ylabel('Shots')
        ax = axes[1, col]
        agg = np.bincount(dec, weights=cts, minlength=2**k)
        top = np.argsort(agg)[::-1][:15]; top = top[agg[top] > 0]
        if opt_d not in top: top = np.r_[top, opt_d]
        cor = ['#1c9275' if j == opt_d else '#2f6db5' if fe_d[j] else '#d0614c' for j in top]
        ax.bar(range(len(top)), agg[top], color=cor)
        ax.axhline(cts.sum()/2**k, color='k', ls=':', lw=1, label='esperado sob amostragem uniforme')
        ax.set_xticks(range(len(top)), [''.join(str((j >> b) & 1) for b in range(k)) for j in top], rotation=90, fontsize=7)
        ax.set_title('Agregado por portfólio x1…xk (todas as folgas)', fontsize=9); ax.set_ylabel('Shots'); ax.legend(fontsize=7)
    fig.suptitle(f"{esp['tipo']}, {esp['k']} controles, p={p}. Verde: portfólio ótimo clássico; azul: factível; vermelho: inviável", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def _serie(runs, ks, f):
    m, s = [], []
    for k in ks:
        v = [f(r) for r in runs if r['k'] == k]
        v = [x for x in v if x is not None]
        m.append(np.mean(v) if v else np.nan); s.append(np.std(v) if v else np.nan)
    return np.array(m), np.array(s)


def figura_escalabilidade(runs, esps, path, titulo):
    ks = sorted({r['k'] for r in runs})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    def plot(ax, f, label, cor, ls='-', marker='o'):
        m, s = _serie(runs, ks, f)
        ok = ~np.isnan(m)
        ax.errorbar(np.array(ks)[ok], m[ok], yerr=s[ok], marker=marker, ls=ls, color=cor, label=label, capsize=3, ms=4)
    ax = axes[0, 0]; plot(ax, lambda r: r['n'], 'qubits lógicos', CORES['E1'])
    plot(ax, lambda r: r['recursos']['fake_kingston']['qubits'], 'qubits físicos (FakeKingston)', CORES['E2'], '--', 's')
    ax.set(ylabel='Qubits')
    ax = axes[0, 1]; plot(ax, lambda r: r['recursos']['logico']['profundidade'], 'lógico', CORES['E1'])
    plot(ax, lambda r: r['recursos']['fake_kingston']['profundidade'], 'transpilado FakeKingston', CORES['E2'], '--', 's')
    ax.set(ylabel='Profundidade do circuito', yscale='log')
    ax = axes[0, 2]; plot(ax, lambda r: r['recursos']['logico']['gates_2q'], 'RZZ lógicas (= #ZZ)', CORES['E1'])
    plot(ax, lambda r: r['recursos']['base_nativa_sem_roteamento']['gates_2q'], 'CZ sem roteamento', CORES['unif'], ':', '^')
    plot(ax, lambda r: r['recursos']['fake_kingston']['gates_2q'], 'CZ FakeKingston', CORES['E2'], '--', 's')
    ax.set(ylabel='Gates de dois qubits', yscale='log')
    for ax, key, lab in ((axes[1, 0], 'Popt', 'Popt'), (axes[1, 1], 'Pfeas', 'Pfeas')):
        plot(ax, lambda r, key=key: r['E1']['exata'][key], 'E1 ideal (exata)', CORES['E1'])
        plot(ax, lambda r, key=key: r['E2']['shots'][key] if 'E2' in r else None, 'E2 FakeKingston', CORES['E2'], '--', 's')
        ax.plot(ks, [esps[k]['uniforme'][key] for k in ks], color=CORES['unif'], ls=':', marker='^', label='amostragem uniforme')
        ax.set(ylabel=lab)
    axes[1, 0].set_yscale('log')
    ax = axes[1, 2]; plot(ax, lambda r: r['otimizacao']['tempo_s'], 'otimização clássica (E1)', CORES['E1'])
    plot(ax, lambda r: r['E2']['tempo_execucao_s'] if 'E2' in r else None, 'execução E2 (simulação ruidosa)', CORES['E2'], '--', 's')
    ax.set(ylabel='Tempo (s)', yscale='log')
    for ax in axes.flat:
        ax.set_xlabel('Tamanho da instância (controles)'); ax.legend(fontsize=7)
        if max(ks) > 7: ax.axvspan(7.5, max(ks)+.5, color='#f2efe6', zorder=0)
    fig.suptitle(titulo + ' (faixa sombreada: extensão sintética; barras: desvio padrão entre sementes)', fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def figura_qualidade(runs_por_cod, esps_por_cod, path):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    estilos = {'folga': ('-', 'o'), 'desbalanceada': ('--', 's')}
    for cod, runs in runs_por_cod.items():
        ks = sorted({r['k'] for r in runs}); ls, mk = estilos[cod]
        for env, cor in (('E1', CORES['E1']), ('E2', CORES['E2'])):
            src = 'exata' if env == 'E1' else 'shots'
            m, s = _serie(runs, ks, lambda r: r[env][src]['ganho95_vs_uniforme'] if env in r else None)
            ok = ~np.isnan(m)
            axes[0].errorbar(np.array(ks)[ok], m[ok], yerr=s[ok], color=cor, ls=ls, marker=mk, capsize=3, ms=4, label=f'{cod} {env}')
            m, s = _serie(runs, ks, lambda r: r[env]['shots']['R_app'] if env in r else None)
            ok = ~np.isnan(m)
            axes[1].errorbar(np.array(ks)[ok], m[ok], yerr=s[ok], color=cor, ls=ls, marker=mk, capsize=3, ms=4, label=f'{cod} {env}')
            m, s = _serie(runs, ks, lambda r: (r[env][src]['Pqubo'] if env in r else None))
            ok = ~np.isnan(m)
            axes[2].errorbar(np.array(ks)[ok], m[ok], yerr=s[ok], color=cor, ls=ls, marker=mk, capsize=3, ms=4, label=f'{cod} {env} (codificação)')
        axes[2].plot(ks, [esps_por_cod[cod][k]['uniforme']['Pqubo'] for k in ks], color=CORES['unif'], ls=ls, label=f'{cod} uniforme')
    axes[0].axhline(1, color='k', lw=.8); axes[0].set(ylabel='P95 QAOA / P95 uniforme', title='Ganho sobre amostragem aleatória')
    axes[1].axhline(LIMIAR, color='k', lw=.8, ls=':'); axes[1].set(ylabel='V(melhor factível amostrado)/V*', title='Qualidade da melhor decisão')
    axes[2].set(ylabel='P(estado fundamental exato do QUBO)', yscale='log', title='Fidelidade da codificação')
    for ax in axes: ax.set_xlabel('Controles'); ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def figura_p(runs, path):
    ps = sorted({r['p'] for r in runs})
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.5))
    def bar(ax, f, title, log=False):
        for j, (env, cor) in enumerate((('E1', CORES['E1']), ('E2', CORES['E2']))):
            vals = [[f(r, env) for r in runs if r['p'] == p] for p in ps]
            if any(v is None for vv in vals for v in vv): continue
            ax.bar(np.array(ps)+(j-.5)*.35, [np.mean(v) for v in vals], .35, yerr=[np.std(v) for v in vals], color=cor, label=env, capsize=3)
        ax.set(title=title, xlabel='p', xticks=ps)
        if log: ax.set_yscale('log')
        ax.legend(fontsize=7)
    src = lambda env: 'exata' if env == 'E1' else 'shots'
    bar(axes[0], lambda r, e: r[e][src(e)]['energia_esperada'] if e in r else None, '⟨H_C⟩ final')
    bar(axes[1], lambda r, e: r[e][src(e)]['Popt'] if e in r else None, 'Popt')
    bar(axes[2], lambda r, e: r[e][src(e)]['Pfeas'] if e in r else None, 'Pfeas')
    bar(axes[3], lambda r, e: r['recursos']['logico' if e == 'E1' else 'fake_kingston']['gates_2q'], 'Gates 2q (E1 lógico, E2 ISA)')
    bar(axes[4], lambda r, e: r['otimizacao']['tempo_s'] if e == 'E1' else r['E2']['tempo_execucao_s'], 'Tempo (s): otimização / execução E2')
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def figuras_circuito(esp, params, med, tc, pasta):
    escala = max(float(np.max(np.abs(esp['E']))), 1.)
    dqc = circuito_didatico(esp['h'], esp['J'], params, escala)
    try:
        fig = dqc.draw('mpl', fold=-1 if esp['n'] <= 5 else 30, scale=.7)
        fig.savefig(pasta/f"circuito_logico_{esp['tipo']}_k{esp['k']}.png", dpi=120, bbox_inches='tight'); plt.close(fig)
        fig = tc.draw('mpl', idle_wires=False, fold=45, scale=.6)
        fig.savefig(pasta/f"circuito_transpilado_{esp['tipo']}_k{esp['k']}.png", dpi=110, bbox_inches='tight'); plt.close(fig)
    except Exception as e:  # pylatexenc ausente: mantém desenho em texto
        print('figura mpl indisponível:', e)
    (pasta/f"circuito_logico_{esp['tipo']}_k{esp['k']}.txt").write_text(str(dqc.draw('text', fold=120)), encoding='utf-8')
    (pasta/f"circuito_transpilado_{esp['tipo']}_k{esp['k']}.txt").write_text(str(tc.draw('text', idle_wires=False, fold=160)), encoding='utf-8')


# ----------------------------------------------------------------------------- estudos complementares
def estudo_ruido(esp, params, med, args):
    """Seção 19: separa transpilação, leitura, porta e decoerência no mesmo circuito ISA."""
    from qiskit_aer.noise import NoiseModel
    tc = transpilar(med, args.nivel, 42)
    be = fake_kingston()
    variantes = {'ISA_sem_ruido': None,
                 'so_leitura': dict(gate_error=False, readout_error=True, thermal_relaxation=False),
                 'so_decoerencia_T1T2': dict(gate_error=False, readout_error=False, thermal_relaxation=True),
                 'so_erro_porta_despolarizante': dict(gate_error=True, readout_error=False, thermal_relaxation=False),
                 'completo_FakeKingston': dict(gate_error=True, readout_error=True, thermal_relaxation=True)}
    out = {}
    for nome, cfg in variantes.items():
        nm = NoiseModel() if cfg is None else NoiseModel.from_backend(be, **cfg)
        t = time.perf_counter(); cts = simular_ruido(tc, args.shots, 42, nm); dt = time.perf_counter()-t
        m = metricas(esp, counts=cts)
        out[nome] = {k: m[k] for k in ('Popt', 'Pfeas', 'P95', 'Pqubo', 'R_app', 'energia_esperada')} | {'tempo_s': dt}
    prob = np.abs(estado(params, esp['E'], esp['n']))**2
    out['E1_ideal_exata'] = {k: float(prob[esp[mk]].sum()) for k, mk in (('Popt', 'opt'), ('Pfeas', 'feas'), ('P95', 'good'), ('Pqubo', 'ground'))}
    return out


def estudo_transpilacao(esp, med, args):
    """Seções 17/20: logical x base nativa x FakeKingston em níveis 0–3."""
    out = {'logico': recursos(med)}
    for nivel in (0, 1, 2, 3):
        out[f'sem_roteamento_nivel{nivel}'] = recursos(transpilar(med, nivel, 42, roteado=False))
        t = time.perf_counter(); tc = transpilar(med, nivel, 42); dt = time.perf_counter()-t
        m = metricas(esp, counts=simular_ruido(tc, args.shots, 42))
        out[f'fake_kingston_nivel{nivel}'] = recursos(tc) | {'tempo_transpilacao_s': dt} | {k: m[k] for k in ('Popt', 'Pfeas', 'P95')}
    return out


# ----------------------------------------------------------------------------- tabelas
def fmt(x, nd=4):
    if x is None: return '–'
    if isinstance(x, (bool, np.bool_)): return 'sim' if x else 'não'
    if isinstance(x, (int, np.integer)): return str(int(x))
    if isinstance(x, float) and (math.isnan(x)): return '–'
    return f'{x:.{nd}f}'


def media(vals):
    vals = [v for v in vals if v is not None]
    return (float(np.mean(vals)), float(np.std(vals))) if vals else (None, None)


def ms(vals, nd=4, pct=False):
    m, s = media(vals)
    if m is None: return '–'
    if pct: return f'{100*m:.{max(nd-2,1)}f}% ± {100*s:.{max(nd-2,1)}f}'
    return f'{m:.{nd}f} ± {s:.{nd}f}'


def salvar_modelo(esp, pasta):
    pasta.mkdir(parents=True, exist_ok=True)
    nomes = [c['nome'] for c in esp['d']['controles']] + [f'folga_2^{m}' for m in range(esp['n']-esp['k'])]
    with (pasta/'ising_h.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['qubit', 'variavel', 'h_i'])
        for i, v in enumerate(esp['h']): w.writerow([i, nomes[i], v])
    with (pasta/'ising_J.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['i', 'j', 'J_ij'])
        for i in range(esp['n']):
            for j in range(i+1, esp['n']):
                if abs(esp['J'][i, j]) > 1e-12: w.writerow([i, j, esp['J'][i, j]])
    np.savetxt(pasta/'qubo_Q.csv', esp['Q'], delimiter=',', fmt='%.10g')
    ordem = np.argsort(esp['E'], kind='stable')[:20]
    with (pasta/'espectro_20_menores.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['posicao', 'bitstring', 'decisao', 'energia', 'factivel', 'otimo', 'valor', 'custo', 'cobertura'])
        for pos, i in enumerate(ordem, 1):
            dd = decodificar(esp, i)
            w.writerow([pos, dd['bitstring'], dd['decisao'], dd['energia'], dd['factivel'], dd['otimo'], dd['valor'], dd['custo'], dd['cobertura']])
    (pasta/'instancia.json').write_text(json.dumps(esp['d'], indent=2, ensure_ascii=False), encoding='utf-8')
    (pasta/'caracteristicas.json').write_text(json.dumps(caracteristicas(esp) | {'C_ising': esp['C'], 'offset_qubo': esp['c0'],
        'erro_max_ising_vs_qubo': esp['erro_ising'], 'erro_max_matriz_vs_energia': esp['erro_qubo'],
        'fundamental_e_otimo_original': esp['ground_e_otimo'], 'uniforme': esp['uniforme']}, indent=2, ensure_ascii=False), encoding='utf-8')
    return [(int(i), decodificar(esp, i)) for i in ordem]


def linhas_csv(runs):
    linhas = []
    for r in runs:
        for env in ('E1', 'E2'):
            if env not in r: continue
            m = r[env]['exata'] if env == 'E1' else r[env]['shots']
            sh = r[env]['shots']
            rec = r['recursos']['logico' if env == 'E1' else 'fake_kingston']
            rp = dict(r[env]['reprodutibilidade']); rp['penalty_lambda'] = json.dumps(rp['penalty_lambda'])
            rp['gamma'] = json.dumps(rp['gamma']); rp['beta'] = json.dumps(rp['beta'])
            linhas.append(rp | {'E_star': r['E_otima'], 'E_QAOA': sh['E_QAOA'], 'energia_esperada': m['energia_esperada'],
                                'Popt': m['Popt'], 'Pfeas': m['Pfeas'], 'Pqubo': m['Pqubo'], 'P95': m['P95'], 'Paceitavel': m['Paceitavel'],
                                'probabilidades': m['fonte'], 'R_app': sh['R_app'], 'ganho95_vs_uniforme': m['ganho95_vs_uniforme'],
                                'util': m['util'] if env == 'E1' else sh['util'],
                                'melhor_bitstring': sh['melhor_energia']['bitstring'], 'melhor_decisao_factivel': (sh['melhor_factivel'] or {}).get('decisao'),
                                'valor_melhor_factivel': (sh['melhor_factivel'] or {}).get('valor'),
                                'menor_energia_e_mais_frequente': sh['menor_energia_e_mais_frequente'],
                                'profundidade': rec['profundidade'], 'gates': rec['gates'], 'gates_2q': rec['gates_2q'],
                                'interacoes_ZZ': r['recursos']['logico']['gates_2q'],
                                'tempo_otimizacao_s': r['otimizacao']['tempo_s'], 'tempo_execucao_s': r[env]['tempo_execucao_s']})
    return linhas


def utilidade(runs, ks, env):
    src = 'exata' if env == 'E1' else 'shots'
    res = {}
    for k in ks:
        rs = [r for r in runs if r['k'] == k and env in r]
        if rs:
            a = sum(r[env]['shots']['decisao_util'] for r in rs); b = sum(r[env][src]['util'] for r in rs)
            res[k] = {'sementes_decisao_util': a, 'sementes_uteis': b, 'sementes': len(rs),
                      'decisao_util': a >= (len(rs)//2+1), 'util': b >= (len(rs)//2+1)}
    return res


def gerar_tabelas(resumo, out):
    L = ['# Tabelas geradas por `exercicio.py`', '',
         f"Configuração: sementes {resumo['config']['seeds']}, {resumo['config']['shots']} shots, COBYLA com até {resumo['config']['maxiter']}·p avaliações por início, "
         f"{resumo['config']['restarts']} inícios, transpilação FakeKingston nível {resumo['config']['nivel']}. Médias ± desvio padrão entre sementes. "
         'E1: probabilidades exatas do statevector. E2: frequências de shots no FakeKingston. Controles 8+ são sintéticos.', '']
    for cod, bloco in resumo['codificacoes'].items():
        L += [f'## Codificação: {cod}', '', '### Características do problema (seção 12.1)', '',
              '| Controles | Instância | Var. decisão | Auxiliares | Qubits | Orçamento B | Termos lineares | Termos quadráticos | #ZZ | Densidade | Portfólios factíveis | max abs(E)/V* | Penalidade |',
              '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
        for k, c in bloco['caracteristicas'].items():
            pen = ', '.join(f'{a}={b:.4g}' for a, b in c.items() if a in ('rho', 'lambda1', 'lambda2'))
            L.append(f"| {k} | {nome_instancia(int(k))} | {c['variaveis_decisao']} | {c['variaveis_auxiliares']} | {c['variaveis_total']} | {c['orcamento']} | {c['termos_lineares']} | {c['termos_quadraticos']} | {c['interacoes_ZZ']} | {c['densidade_qubo']:.2f} | {c['portfolios_factiveis']}/{c['portfolios_totais']} | {c['razao_escala_valor_otimo']:.0f} | {pen} |")
        L += ['', '### Tabela 2 — Comparação entre solução clássica e QAOA (p=1)', '',
              '| Controles | Qubits | Clássico (decisão) | QAOA E1 melhor factível | QAOA E2 melhor factível | E* | E_QAOA E1 | E_QAOA E2 | Popt E1 | Popt E2 | Pfeas E1 | Pfeas E2 | Popt uniforme | Pfeas uniforme |',
              '|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for k, t in bloco['tabela2'].items():
            L.append(f"| {k} | {t['qubits']} | `{t['classico']}` | {t['qaoa_E1']} | {t['qaoa_E2']} | {t['E_star']:.4f} | {t['E_QAOA_E1']} | {t['E_QAOA_E2']} | {t['Popt_E1']} | {t['Popt_E2']} | {t['Pfeas_E1']} | {t['Pfeas_E2']} | {t['Popt_unif']:.4f} | {t['Pfeas_unif']:.4f} |")
        L += ['', 'Decisões no formato x_1…x_k (1 = controle selecionado); entre parênteses, quantas sementes obtiveram aquela decisão como melhor amostra factível. E_QAOA é a menor energia QUBO amostrada.', '',
              '### Recursos quânticos e utilidade (p=1)', '',
              '| Controles | Qubits | Depth lógico | Depth ISA | Gates lógicos | Gates ISA | 2q lógicos (#ZZ) | CZ sem roteamento | CZ ISA | Shots | P95 E1/unif | P95 E2/unif | R_app E1 | R_app E2 | Útil E1 (A, B) | Útil E2 (A, B) | t otim. (s) | t E2 (s) |',
              '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|']
        for k, t in bloco['recursos'].items():
            L.append('| ' + ' | '.join(str(t[c]) for c in ['k', 'qubits', 'depth_log', 'depth_isa', 'gates_log', 'gates_isa', 'zz', 'cz_base', 'cz_isa', 'shots', 'g95_E1', 'g95_E2', 'rapp_E1', 'rapp_E2', 'util_E1', 'util_E2', 't_opt', 't_E2']) + ' |')
        if bloco.get('profundidades'):
            L += ['', f"### Influência de p ({bloco['profundidades']['k']} controles)", '',
                  '| p | Ambiente | ⟨H_C⟩ | Popt | Pfeas | P95 | Pqubo | Depth | Gates 2q | t otimização (s) | Avaliações |', '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            for row in bloco['profundidades']['linhas']:
                L.append('| ' + ' | '.join(str(x) for x in row) + ' |')
        if bloco.get('ruido'):
            L += ['', '### Efeito do ruído por tamanho (p=1)', '',
                  '| Controles | Popt ideal | Popt ruído | Pfeas ideal | Pfeas ruído | P95 ideal | P95 ruído | Δ_ruído Popt | Δ_ruído Pfeas | CZ ISA |', '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
            for row in bloco['ruido']:
                L.append('| ' + ' | '.join(str(x) for x in row) + ' |')
        L += ['', '### Benefício no contexto do ODS (maior instância da base real e maior com E2)', '']
        for k, cen in bloco['beneficio'].items():
            L += [f'**{k} controles** (referência = implantar na ordem numérica CIS até esgotar o orçamento)', '',
                  '| Cenário | Valor de segurança | Δ valor vs ref. | Cobertura ATT&CK | Δ cobertura (p.p.) | Custo | % do ótimo |', '|---|---:|---:|---:|---:|---:|---:|']
            for nome, c in cen.items():
                L.append(f"| {nome} | {c['valor']:.4f} | {c['delta_valor_pct_vs_ref']:+.2f}% | {100*c['cobertura']:.1f}% | {c['delta_cobertura_pp_vs_ref']:+.1f} | {c['custo']:.2f} | {100*c['fracao_do_otimo']:.1f}% |")
            L.append('')
    if resumo.get('estudo_ruido'):
        L += ['## Decomposição do ruído (7 controles, folga, p=1, semente 42)', '', '| Variante | Popt | Pfeas | P95 | Pqubo | R_app | ⟨H_C⟩ |', '|---|---:|---:|---:|---:|---:|---:|']
        for nome, m in resumo['estudo_ruido'].items():
            L.append(f"| {nome} | {m['Popt']:.4f} | {m['Pfeas']:.4f} | {m['P95']:.4f} | {m['Pqubo']:.4f} | {fmt(m.get('R_app'))} | {fmt(m.get('energia_esperada'), 1)} |")
        L.append('')
    if resumo.get('estudo_transpilacao'):
        L += ['## Circuito antes e depois da transpilação (7 controles, folga, p=1)', '', '| Versão | Qubits | Depth | Gates | Two-qubit gates | Popt E2 | Pfeas E2 |', '|---|---:|---:|---:|---:|---:|---:|']
        for nome, m in resumo['estudo_transpilacao'].items():
            L.append(f"| {nome} | {m['qubits']} | {m['profundidade']} | {m['gates']} | {m['gates_2q']} | {fmt(m.get('Popt'))} | {fmt(m.get('Pfeas'))} |")
        L.append('')
    (out/'TABELAS.md').write_text('\n'.join(L), encoding='utf-8')


def consolidar(cod, runs_cod, esps, ks, kp):
    p1 = [r for r in runs_cod if r['p'] == 1]
    bloco = {'caracteristicas': {k: caracteristicas(esps[k]) for k in ks}, 'tabela2': {}, 'recursos': {}, 'ruido': [], 'beneficio': {},
             'utilidade': {'E1': utilidade(p1, ks, 'E1'), 'E2': utilidade(p1, ks, 'E2')}}
    for k in ks:
        rs = [r for r in p1 if r['k'] == k]; e2 = [r for r in rs if 'E2' in r]
        def decisoes(env):
            if env == 'E2' and not e2: return '–'
            dd = {}
            for r in (rs if env == 'E1' else e2):
                b = r[env]['shots']['melhor_factivel']
                key = b['decisao'] if b else 'nenhuma'
                dd[key] = dd.get(key, 0)+1
            ok = esps[k]['bits'][np.flatnonzero(esps[k]['opt'])[0], :k].astype(int)
            ok = ''.join(map(str, ok))
            return ', '.join(f"`{x}`{' ✓' if x == ok else ''} ({c})" for x, c in dd.items())
        bloco['tabela2'][k] = {'qubits': esps[k]['n'], 'classico': ''.join(map(str, esps[k]['bits'][np.flatnonzero(esps[k]['opt'])[0], :k].astype(int))),
                               'qaoa_E1': decisoes('E1'), 'qaoa_E2': decisoes('E2'), 'E_star': float(esps[k]['E'].min()),
                               'E_QAOA_E1': ms([r['E1']['shots']['E_QAOA'] for r in rs]), 'E_QAOA_E2': ms([r['E2']['shots']['E_QAOA'] for r in e2]) if e2 else '–',
                               'Popt_E1': ms([r['E1']['exata']['Popt'] for r in rs]), 'Popt_E2': ms([r['E2']['shots']['Popt'] for r in e2]) if e2 else '–',
                               'Pfeas_E1': ms([r['E1']['exata']['Pfeas'] for r in rs]), 'Pfeas_E2': ms([r['E2']['shots']['Pfeas'] for r in e2]) if e2 else '–',
                               'Popt_unif': esps[k]['uniforme']['Popt'], 'Pfeas_unif': esps[k]['uniforme']['Pfeas']}
        rec = lambda tipo, key, sub=rs: media([r['recursos'][tipo][key] for r in sub])[0]
        u1, u2 = bloco['utilidade']['E1'].get(k), bloco['utilidade']['E2'].get(k)
        bloco['recursos'][k] = {'k': k, 'qubits': esps[k]['n'], 'depth_log': f"{rec('logico', 'profundidade'):.0f}", 'depth_isa': f"{rec('fake_kingston', 'profundidade'):.0f}",
                                'gates_log': f"{rec('logico', 'gates'):.0f}", 'gates_isa': f"{rec('fake_kingston', 'gates'):.0f}", 'zz': f"{rec('logico', 'gates_2q'):.0f}",
                                'cz_base': f"{rec('base_nativa_sem_roteamento', 'gates_2q'):.0f}", 'cz_isa': f"{rec('fake_kingston', 'gates_2q'):.0f}",
                                'shots': rs[0]['E1']['shots']['shots'],
                                'g95_E1': ms([r['E1']['exata']['ganho95_vs_uniforme'] for r in rs], 2),
                                'g95_E2': ms([r['E2']['shots']['ganho95_vs_uniforme'] for r in e2], 2) if e2 else '–',
                                'rapp_E1': ms([r['E1']['shots']['R_app'] for r in rs], 3), 'rapp_E2': ms([r['E2']['shots']['R_app'] for r in e2], 3) if e2 else '–',
                                'util_E1': f"A {u1['sementes_decisao_util']}/{u1['sementes']}, B {u1['sementes_uteis']}/{u1['sementes']}" if u1 else '–',
                                'util_E2': f"A {u2['sementes_decisao_util']}/{u2['sementes']}, B {u2['sementes_uteis']}/{u2['sementes']}" if u2 else 'não executado',
                                't_opt': f"{media([r['otimizacao']['tempo_s'] for r in rs])[0]:.2f}",
                                't_E2': f"{media([r['E2']['tempo_execucao_s'] for r in e2])[0]:.1f}" if e2 else '–'}
        if e2:
            bloco['ruido'].append([k, ms([r['E1']['exata']['Popt'] for r in e2]), ms([r['E2']['shots']['Popt'] for r in e2]),
                                   ms([r['E1']['exata']['Pfeas'] for r in e2]), ms([r['E2']['shots']['Pfeas'] for r in e2]),
                                   ms([r['E1']['exata']['P95'] for r in e2]), ms([r['E2']['shots']['P95'] for r in e2]),
                                   ms([r['delta_ruido']['Popt'] for r in e2]), ms([r['delta_ruido']['Pfeas'] for r in e2]),
                                   f"{media([r['recursos']['fake_kingston']['gates_2q'] for r in e2])[0]:.0f}"])
    alvo = [k for k in ks if k <= 7]
    alvo = {max(alvo)} if alvo else set()
    com_e2 = [k for k in ks if any('E2' in r for r in p1 if r['k'] == k)]
    if com_e2: alvo.add(max(com_e2))
    for k in sorted(alvo):
        bloco['beneficio'][k] = beneficio(esps[k], [r for r in p1 if r['k'] == k])
    pr = [r for r in runs_cod if r['k'] == kp]
    if len({r['p'] for r in pr}) > 1:
        linhas = []
        for p in sorted({r['p'] for r in pr}):
            g = [r for r in pr if r['p'] == p]
            for env in ('E1', 'E2'):
                gg = [r for r in g if env in r]
                if not gg: continue
                src = 'exata' if env == 'E1' else 'shots'
                rk = 'logico' if env == 'E1' else 'fake_kingston'
                linhas.append([p, env, ms([r[env][src]['energia_esperada'] for r in gg], 1), ms([r[env][src]['Popt'] for r in gg]),
                               ms([r[env][src]['Pfeas'] for r in gg]), ms([r[env][src]['P95'] for r in gg]), ms([r[env][src]['Pqubo'] for r in gg]),
                               f"{media([r['recursos'][rk]['profundidade'] for r in gg])[0]:.0f}", f"{media([r['recursos'][rk]['gates_2q'] for r in gg])[0]:.0f}",
                               f"{media([r['otimizacao']['tempo_s'] for r in gg])[0]:.2f}", f"{media([r['otimizacao']['avaliacoes_total'] for r in gg])[0]:.0f}"])
        bloco['profundidades'] = {'k': kp, 'linhas': linhas}
    bloco['maior_instancia_util'] = {}
    for env in ('E1', 'E2'):
        for nivel, chave in (('A_decisao_util', 'decisao_util'), ('B_atribuivel_QAOA', 'util')):
            v = [k for k, u in bloco['utilidade'][env].items() if u[chave]]
            bloco['maior_instancia_util'][f'{env}_{nivel}'] = max(v) if v else None
    return bloco


# ----------------------------------------------------------------------------- principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dados', type=Path, default=ROOT/'instancia.json')
    ap.add_argument('--saida', type=Path, default=ROOT/'resultados_exercicio')
    ap.add_argument('--tamanhos', type=int, nargs='+', default=[3, 4, 5, 6, 7, 9, 11, 13, 15])
    ap.add_argument('--tamanho-p', type=int, default=7, help='instância usada para p=1,2,3')
    ap.add_argument('--profundidades', type=int, nargs='+', default=[1, 2, 3])
    ap.add_argument('--seeds', type=int, nargs='+', default=[42, 43, 44])
    ap.add_argument('--shots', type=int, default=2048)
    ap.add_argument('--maxiter', type=int, default=150, help='avaliações COBYLA por início, multiplicado por p')
    ap.add_argument('--restarts', type=int, default=2)
    ap.add_argument('--nivel', type=int, default=2, help='optimization_level da transpilação FakeKingston')
    ap.add_argument('--max-qubits-e2', type=int, default=17, help='E2 só até este número de qubits (tempo de simulação ruidosa)')
    ap.add_argument('--codificacoes', nargs='+', default=['folga', 'desbalanceada'], choices=['folga', 'desbalanceada'])
    ap.add_argument('--calibracao', type=int, nargs='+', default=[3, 4, 5], help='instâncias que calibram (l1,l2)')
    ap.add_argument('--sem-estudos', action='store_true', help='pula decomposição do ruído e níveis de transpilação')
    args = ap.parse_args()
    import qiskit, qiskit_aer, qiskit_ibm_runtime, scipy
    base = json.loads(args.dados.read_text())
    controles = familia(base, max(args.tamanhos))
    args.saida.mkdir(parents=True, exist_ok=True)
    resumo = {'config': {**{k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                         'versoes': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                                     'qiskit': qiskit.__version__, 'qiskit_aer': qiskit_aer.__version__,
                                     'qiskit_ibm_runtime': qiskit_ibm_runtime.__version__},
                         'backend_E2': {'nome': fake_kingston().name, 'qubits': fake_kingston().num_qubits,
                                        'portas': sorted(fake_kingston().operation_names)},
                         'limiar_qualidade': LIMIAR, 'ods': ODS, 'problema': PROBLEMA,
                         'criterio_util': {'A_decisao_util': 'R_app >= 0.95: o melhor portfólio factível amostrado atinge 95% do ótimo',
                                           'B_atribuivel_QAOA': 'A e P95 do QAOA (exata em E1; limite inferior de Wilson 95% em E2) > P95 da amostragem uniforme',
                                           'instancia': 'útil no nível se a maioria das sementes satisfaz'}},
              'codificacoes': {}}
    (args.saida/'controles_familia.json').write_text(json.dumps(controles, indent=2, ensure_ascii=False), encoding='utf-8')
    todas = []
    for cod in args.codificacoes:
        lam = None
        if cod == 'desbalanceada':
            lam, cal = calibrar_desbalanceada(base, controles, args.calibracao)
            resumo['calibracao_desbalanceada'] = cal
            print('calibração desbalanceada:', cal['escolhido'], flush=True)
        esps, runs = {}, []
        for k in args.tamanhos:
            esp = codificar(instancia(base, controles, k), cod, lam); esps[k] = esp
            salvar_modelo(esp, args.saida/'modelos'/f'{cod}_k{k}')
            ps = args.profundidades if k == args.tamanho_p else [1]
            e2 = esp['n'] <= args.max_qubits_e2
            for seed in args.seeds:
                prev = None
                for p in ps:
                    out = args.saida/'execucoes'/f'{cod}_k{k}_p{p}_s{seed}'
                    r, params, tc, med = executar(esp, p, seed, args, interp(prev[:p-1]).tolist() + interp(prev[p-1:]).tolist() if prev is not None else None, out, e2)
                    prev = params; runs.append(r)
                    if seed == args.seeds[0] and p == 1 and k in (min(args.tamanhos), args.tamanho_p):
                        figuras_circuito(esp, params, med, tc, args.saida)
                    if (not args.sem_estudos and cod == 'folga' and k == args.tamanho_p and p == 1 and seed == 42):
                        resumo['estudo_ruido'] = estudo_ruido(esp, params, med, args)
                        resumo['estudo_transpilacao'] = estudo_transpilacao(esp, med, args)
                    e2s = r['E2']['shots'] if 'E2' in r else None
                    print(f"{cod} k={k} n={esp['n']} p={p} s={seed}: Popt E1={r['E1']['exata']['Popt']:.4f} "
                          f"P95/unif={r['E1']['exata']['ganho95_vs_uniforme']:.2f}"
                          + (f" | E2 Popt={e2s['Popt']:.4f} Pfeas={e2s['Pfeas']:.3f}" if e2s else '')
                          + f" | {r['tempo_total_s']:.1f}s", flush=True)
        resumo['codificacoes'][cod] = consolidar(cod, runs, esps, args.tamanhos, args.tamanho_p)
        figura_escalabilidade([r for r in runs if r['p'] == 1], esps, args.saida/f'escalabilidade_{cod}.png',
                              f'Escalabilidade, codificação {cod}, p=1')
        if len(args.profundidades) > 1:
            figura_p([r for r in runs if r['k'] == args.tamanho_p], args.saida/f'profundidade_p_{cod}.png')
        todas.append((cod, runs, esps))
    figura_qualidade({c: [r for r in rs if r['p'] == 1] for c, rs, _ in todas}, {c: e for c, _, e in todas}, args.saida/'qualidade_e_utilidade.png')
    linhas = [l for _, rs, _ in todas for l in linhas_csv(rs)]
    with (args.saida/'execucoes.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0])); w.writeheader(); w.writerows(linhas)
    (args.saida/'execucoes.json').write_text(json.dumps([r for _, rs, _ in todas for r in rs], indent=1, ensure_ascii=False), encoding='utf-8')
    (args.saida/'resumo.json').write_text(json.dumps(resumo, indent=2, ensure_ascii=False, default=str), encoding='utf-8')
    gerar_tabelas(resumo, args.saida)
    print('maior instância útil:', {c: b['maior_instancia_util'] for c, b in resumo['codificacoes'].items()})


if __name__ == '__main__':
    main()
