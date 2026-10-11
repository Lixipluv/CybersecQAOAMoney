"""Seleção de controles: modelo do artigo, QUBO, Ising e QAOA."""
import argparse
import csv
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def validar(d):
    n = len(d['controles'])
    if n == 0 or not isinstance(d['orcamento'], int) or d['orcamento'] < 0:
        raise ValueError('Orçamento deve ser inteiro não negativo; forneça controles.')
    for c in d['controles']:
        if not isinstance(c['custo'], int) or c['custo'] <= 0:
            raise ValueError('Custos devem ser inteiros positivos na mesma unidade.')
        for k in ['eficacia', 'adesao', 'prioridade', 'impacto']:
            if not 0 <= c[k] <= 1:
                raise ValueError(f'{k} fora de [0,1]')
        if len(c['ativos']) != len(d['ativos']) or len(c['taticas']) != len(d['relevancia_taticas']):
            raise ValueError('Dimensão de ativos/táticas incompatível.')
        if any(v not in (0, 1) for v in c['ativos'] + c['taticas']):
            raise ValueError('Relações devem ser binárias.')
    for key in ['pesos_valor', 'pesos_impactos']:
        if any(v < 0 for v in d[key]) or not np.isclose(sum(d[key]), 1):
            raise ValueError(f'{key} deve somar 1 com pesos não negativos.')
    if len(d['pesos_valor']) != 4 or any(len(a['impactos']) != len(d['pesos_impactos']) for a in d['ativos']):
        raise ValueError('Dimensão dos pesos incorreta.')
    for a in d['ativos']:
        if any(not 0 <= v <= 1 for v in a['impactos'] + [a['exposicao']]):
            raise ValueError('Impactos/exposição fora de [0,1].')
    dep = np.array(d['dependencias'], float)
    if dep.shape != (len(d['ativos']), len(d['ativos'])) or np.any((dep < 0) | (dep > 1)):
        raise ValueError('Matriz de dependências inválida.')
    if any(not 0 <= d[k] <= 1 for k in ['eta', 'lambda_prioridade', 'cobertura_minima']):
        raise ValueError('Parâmetro fora de [0,1].')
    q = np.array(d['relevancia_taticas'])
    if np.any((q < 0) | (q > 1)) or q.sum() <= 0:
        raise ValueError('Relevâncias inválidas.')


def valores(d):
    validar(d)
    v = np.array([a['impactos'] for a in d['ativos']]) @ d['pesos_impactos']
    dep = np.array(d['dependencias'], float)
    np.fill_diagonal(dep, 0)
    v = (v + d['eta'] * dep @ v) / (1 + d['eta'] * dep.sum(axis=1))
    r = np.array([a['exposicao'] for a in d['ativos']])
    q = np.array(d['relevancia_taticas'])
    result = []
    for c in d['controles']:
        z = np.array(c['ativos'])
        u = np.dot(z, v*r) / z.sum() if z.sum() else 0
        a = np.dot(q, c['taticas']) / q.sum()
        w = d['lambda_prioridade']*c['prioridade'] + (1-d['lambda_prioridade'])*u
        result.append(np.dot(d['pesos_valor'], [c['eficacia']*c['adesao'], a, w, 1-c['impacto']]))
    return np.array(result)


def modelo(d, rho=None):
    R = valores(d)
    B = d['orcamento']
    m = B.bit_length()
    a = np.array([c['custo'] for c in d['controles']] + [2**j for j in range(m)], float)
    rho = float(R.sum()+1 if rho is None else rho)
    if rho <= R.sum():
        raise ValueError('Use rho > soma(R) para a garantia de factibilidade.')
    linear = rho*(a*a-2*B*a)
    linear[:len(R)] -= R
    Q = np.diag(linear)
    for i in range(len(a)):
        for j in range(i+1, len(a)):
            Q[i,j] = 2*rho*a[i]*a[j]
    # Convenção triangular: b^T Q b + constante.
    h = -linear/2
    J = np.triu(Q, 1)/4
    h -= (J.sum(axis=0)+J.sum(axis=1))
    constant = rho*B*B + linear.sum()/2 + J.sum()
    return R, a, rho, Q, h, J, constant


def tabela(d, rho=None):
    R,a,rho,Q,h,J,C = modelo(d,rho)
    n = len(a)
    if n > 20:
        raise ValueError('Limite de segurança: 20 qubits no simulador exaustivo.')
    indices = np.arange(2**n)
    bits = ((indices[:,None] >> np.arange(n)) & 1).astype(float)
    costs = bits[:,:len(R)] @ a[:len(R)]
    benefits = bits[:,:len(R)] @ R
    E = -benefits + rho*(bits @ a-d['orcamento'])**2
    return R,a,rho,Q,h,J,C,bits,costs,benefits,E


def estado(params, E, n):
    # Energia normalizada melhora condicionamento; gamma refere-se a HC/escala.
    scale = max(float(np.max(np.abs(E))), 1.)
    psi = np.ones(len(E), complex)/np.sqrt(len(E))
    p = len(params)//2
    for gamma,beta in zip(params[:p], params[p:]):
        psi *= np.exp(-1j*gamma*E/scale)
        for q in range(n):
            v = psi.reshape(-1,2,2**q)
            lo,hi = v[:,0,:].copy(), v[:,1,:].copy()
            v[:,0,:] = np.cos(beta)*lo-1j*np.sin(beta)*hi
            v[:,1,:] = np.cos(beta)*hi-1j*np.sin(beta)*lo
    return psi


def circuito(h,J,C,params,E):
    from qiskit import QuantumCircuit
    n,p = len(h),len(params)//2
    scale = max(float(np.max(np.abs(E))),1.)
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for gamma,beta in zip(params[:p],params[p:]):
        for i in range(n):
            if abs(h[i]) > 1e-12: qc.rz(2*gamma*h[i]/scale,i)
        for i in range(n):
            for j in range(i+1,n):
                if abs(J[i,j]) > 1e-12: qc.rzz(2*gamma*J[i,j]/scale,i,j)
        qc.rx(2*beta,range(n))
    qc.global_phase = -sum(params[:p])*C/scale
    return qc


def recursos(qc):
    gates = [i for i in qc.data if i.operation.name not in ['measure','barrier']]
    return {'qubits':qc.num_qubits,'profundidade':qc.depth(), 'gates':len(gates),
            'gates_2q':sum(len(i.qubits)==2 for i in gates)}


def executar(d,p,seed,shots,maxiter,restarts,out,ruido=False,rho=None):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    t0 = time.perf_counter()
    R,a,rho,Q,h,J,C,bits,costs,V,E = tabela(d,rho)
    n,k = len(a),len(R)
    feasible = costs <= d['orcamento']
    best_value = V[feasible].max()
    optimal = feasible & np.isclose(V,best_value,atol=1e-10,rtol=0)
    rng = np.random.default_rng(seed)
    history,runs = [],[]
    opt_start = time.perf_counter()
    for restart in range(restarts):
        initial = np.r_[rng.uniform(0,2*np.pi,p),rng.uniform(0,np.pi,p)]
        trace=[]
        def objective(x):
            energy = float(np.dot(abs(estado(x,E,n))**2,E))
            trace.append(energy)
            return energy/max(float(np.max(np.abs(E))),1.)
        res = minimize(objective, initial, method='COBYLA', options={'maxiter':maxiter,'rhobeg':1.0,'tol':1e-5})
        runs.append({'iniciais':initial.tolist(),'finais':res.x.tolist(),'avaliacoes':int(res.nfev),
                     'energia_inicial':trace[0], 'energia_final':float(res.fun)*max(float(np.max(np.abs(E))),1.),
                     'convergiu':bool(res.success),'mensagem':str(res.message)})
        history.append(trace)
    chosen = min(range(restarts),key=lambda i:runs[i]['energia_final'])
    params = np.array(runs[chosen]['finais'])
    opt_time = time.perf_counter()-opt_start
    prob = abs(estado(params,E,n))**2
    qc = circuito(h,J,C,params,E)
    measured = qc.copy(); measured.measure_all()
    tc = transpile(measured,basis_gates=['rz','sx','x','cx'], coupling_map=CouplingMap.from_line(n),
                   optimization_level=1,seed_transpiler=seed)
    counts = rng.multinomial(shots,prob)
    environment = 'ideal'
    if ruido:
        noise = NoiseModel()
        noise.add_all_qubit_quantum_error(depolarizing_error(.001,1),['sx','x'])
        noise.add_all_qubit_quantum_error(depolarizing_error(.01,2),['cx'])
        noise.add_all_qubit_readout_error(ReadoutError([[.98,.02],[.02,.98]]))
        raw = AerSimulator(noise_model=noise,method='statevector').run(tc,shots=shots,seed_simulator=seed).result().get_counts()
        counts = np.zeros(len(E),int)
        for b,count in raw.items(): counts[int(b.replace(' ',''),2)] += count
        environment = 'ruido_sintetico'
    sampled = np.flatnonzero(counts)
    best = int(sampled[np.argmin(E[sampled])])
    feas_sample = sampled[feasible[sampled]]
    best_feas = int(feas_sample[np.argmax(V[feas_sample])]) if len(feas_sample) else None
    q = np.array(d['relevancia_taticas'])
    coverage = (np.minimum(1,bits[:,:k] @ np.array([c['taticas'] for c in d['controles']])) @ q)/q.sum()
    acceptable = feasible & (V>=d['valor_minimo']) & (coverage>=d['cobertura_minima'])
    mode=int(np.argmax(counts))
    def decode(index):
        if index is None: return None
        return {'bitstring':format(index,f'0{n}b'),'controles':[d['controles'][i]['nome'] for i in range(k) if bits[index,i]],
                'custo':float(costs[index]),'valor':float(V[index]),'cobertura':float(coverage[index]),
                'factivel':bool(feasible[index]),'aceitavel':bool(acceptable[index]),'energia':float(E[index])}
    result = {'controles':k,'auxiliares':n-k,'qubits':n,'restricoes_qubo':1,'p':p,'seed':seed,'shots':shots,
              'ambiente':environment,'rho':rho,'escala_hamiltoniano':max(float(np.max(np.abs(E))),1.),
              'E_otima':float(E.min()),'valor_otimo':float(best_value),'E_melhor_medido':float(E[best]),
              'Popt':float(counts[optimal].sum()/shots),'Pfeas':float(counts[feasible].sum()/shots),
              'Paceitavel':float(counts[acceptable].sum()/shots),
              'Popt_exata_ideal':float(prob[optimal].sum()),'Pfeas_exata_ideal':float(prob[feasible].sum()),
              'energia_esperada_ideal':float(prob@E),'otimizacao_s':opt_time,
              'classico':decode(int(np.flatnonzero(optimal)[0])), 'melhor_medido':decode(best),
              'melhor_factivel':decode(best_feas),'mais_frequente':decode(mode),
              'menor_energia_e_mais_frequente':bool(np.isclose(E[best],E[mode])),
              'aceitabilidade_possivel':bool(acceptable.any()),
              'termos_lineares':int(np.count_nonzero(h)),'interacoes_ZZ':int(np.count_nonzero(J)),
              'logico':recursos(qc),'transpilado':recursos(tc),'otimizacoes':runs}
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'modelo.npz',Q=Q,h=h,J=J,constante=C,offset_qubo=rho*d['orcamento']**2,valores=R)
    (out/'circuito_logico.txt').write_text(str(measured.draw(output='text')),encoding='utf-8')
    (out/'circuito_transpilado.txt').write_text(str(tc.draw(output='text')),encoding='utf-8')
    with (out/'medicoes.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['bitstring','shots','energia','factivel','otimo','valor','custo','cobertura'])
        for i in sampled: w.writerow([format(i,f'0{n}b'),counts[i],E[i],bool(feasible[i]),bool(optimal[i]),V[i],costs[i],coverage[i]])
    fig,ax=plt.subplots(figsize=(8,4))
    for i,trace in enumerate(history): ax.plot(trace,label=f'início {i+1}')
    ax.axhline(E.min(),color='black',ls='--',label='ótimo QUBO')
    ax.set(xlabel='Avaliação da função',ylabel='Energia',title=f'{k} controles, p={p}');ax.legend();fig.tight_layout()
    fig.savefig(out/'energia.png');plt.close(fig)
    top=sampled[np.argsort(counts[sampled])[-15:]]
    fig,ax=plt.subplots(figsize=(10,4));ax.bar([format(i,f'0{n}b') for i in top],counts[top],color=['green' if optimal[i] else 'steelblue' if feasible[i] else 'salmon' for i in top])
    ax.tick_params(axis='x',rotation=90);ax.set(ylabel='Shots',title='Verde: ótimo; azul: factível; vermelho: inviável');fig.tight_layout();fig.savefig(out/'histograma.png');plt.close(fig)
    result['execucao_s']=time.perf_counter()-t0
    (out/'resultado.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dados',type=Path,default=Path(__file__).with_name('instancia.json'))
    parser.add_argument('--saida',type=Path,default=Path(__file__).with_name('resultados'))
    parser.add_argument('--tamanhos',type=int,nargs='+',default=[3,4,5,6,7])
    parser.add_argument('--profundidades',type=int,nargs='+',default=[1,2,3])
    parser.add_argument('--repeticoes',type=int,default=3)
    parser.add_argument('--restarts',type=int,default=2)
    parser.add_argument('--maxiter',type=int,default=150)
    parser.add_argument('--shots',type=int,default=2048)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--ruido',action='store_true')
    args=parser.parse_args()
    if min(args.shots,args.maxiter,args.restarts,args.repeticoes,*args.profundidades,*args.tamanhos)<=0:
        parser.error('Parâmetros devem ser positivos.')
    d=json.loads(args.dados.read_text())
    validar(d)
    if max(args.tamanhos)>len(d['controles']): parser.error('Tamanho excede controles disponíveis.')
    results=[]
    for k in args.tamanhos:
        inst=dict(d);inst['controles']=d['controles'][:k]
        # Orçamento proporcional permite comparar crescimento sem tornar tudo selecionável.
        inst['orcamento']=max(1,int(sum(c['custo'] for c in inst['controles'])*d['fracao_orcamento']))
        ps=args.profundidades if k==max(args.tamanhos) else [args.profundidades[0]]
        for p in ps:
            for rep in range(args.repeticoes):
                for noise in ([False,True] if args.ruido and k==min(args.tamanhos) else [False]):
                    out=args.saida/f'n{k}_p{p}_rep{rep}_{"ruido" if noise else "ideal"}'
                    result=executar(inst,p,args.seed+rep,args.shots,args.maxiter,args.restarts,out,noise)
                    results.append(result)
                    print(f'n={k}, p={p}, {result["ambiente"]}: Popt={result["Popt"]:.3f}, Pfeas={result["Pfeas"]:.3f}',flush=True)
    args.saida.mkdir(parents=True,exist_ok=True)
    (args.saida/'experimentos.json').write_text(json.dumps(results,indent=2,ensure_ascii=False))
    (args.saida/'configuracao.json').write_text(json.dumps({'argumentos':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},'dados':d},indent=2,ensure_ascii=False))
    rows=[{'controles':r['controles'],'p':r['p'],'ambiente':r['ambiente'],'seed':r['seed'],
           **{key:r[key] for key in ['qubits','Popt','Pfeas','Paceitavel','E_otima','E_melhor_medido','otimizacao_s','execucao_s']},
           **{f'{a}_{b}':r[a][b] for a in ['logico','transpilado'] for b in ['profundidade','gates','gates_2q']}} for r in results]
    with (args.saida/'resumo.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    base=[r for r in results if r['p']==args.profundidades[0] and r['ambiente']=='ideal']
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for ax,key,label in zip(axes.flat,['qubits','profundidade','gates_2q','Popt','Pfeas','otimizacao_s'],['Qubits','Profundidade transpilada','Gates 2q transpilados','Popt','Pfeas','Tempo de otimização (s)']):
        for k in args.tamanhos:
            vals=[r['transpilado'][key] if key in ['profundidade','gates_2q'] else r[key] for r in base if r['controles']==k]
            ax.errorbar(k,np.mean(vals),yerr=np.std(vals),fmt='o',color='steelblue')
        ax.set(xlabel='Número de controles',ylabel=label)
    fig.tight_layout();fig.savefig(args.saida/'escalabilidade.png');plt.close(fig)

if __name__=='__main__': main()
