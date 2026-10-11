"""Ablation reproduzível: penalidade, CVaR, escala e warm start heurístico.

A referência exata só participa da avaliação; nunca da seleção de parâmetros.
"""
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from projeto import tabela, valores, estado, circuito, recursos, plt

METHODS = ('baseline', 'penalidade_menor', 'cvar', 'cvar_escala', 'warm_cvar')


def cvar(prob, energies, alpha, order=None):
    """Média ponderada da cauda inferior; inclui massa fracionária na fronteira."""
    if not 0 < alpha <= 1:
        raise ValueError('alpha deve pertencer a (0,1].')
    order = np.argsort(energies, kind='stable') if order is None else order
    p = prob[order]
    mass = np.minimum(p, np.maximum(alpha - (np.cumsum(p)-p), 0))
    return float(mass @ energies[order] / alpha)


def greedy(d, R):
    """Valor/custo, sem acesso à solução exata; empates pelo índice."""
    costs = np.array([c['custo'] for c in d['controles']])
    bits = np.zeros(len(R), int)
    used = 0
    for i in sorted(range(len(R)), key=lambda i: (-R[i]/costs[i], i)):
        if used + costs[i] <= d['orcamento']:
            bits[i] = 1
            used += costs[i]
    slack = d['orcamento'] - used
    return np.r_[bits, [(slack >> m) & 1 for m in range(d['orcamento'].bit_length())]]


def setup(d, method, epsilon=.25):
    R = valores(d)
    # Derivação em ANALISE_MELHORIAS.md; preserva o ótimo para este modelo.
    rho = None if method == 'baseline' else float(R.max()+.05)
    model = tabela(d, rho)
    R,a,rho,Q,h,J,C,bits,costs,V,E = model
    scale = max(float(np.max(np.abs(E))), 1.)
    if method in ('cvar_escala', 'warm_cvar'):
        scale = max(float(np.max(np.abs(h))), float(np.max(np.abs(J))), 1.)
    warm = None
    if method == 'warm_cvar':
        warm = np.clip(greedy(d,R).astype(float),epsilon,1-epsilon)
    return model,scale,warm


def evaluate(prob, model, d, shots, seed):
    R,a,rho,Q,h,J,C,bits,costs,V,E = model
    feasible = costs <= d['orcamento']
    opt_value = float(V[feasible].max())
    optimum = feasible & np.isclose(V,opt_value,rtol=0,atol=1e-10)
    ground = np.isclose(E,E.min(),rtol=0,atol=1e-8)
    q = np.array(d['relevancia_taticas'])
    cov = np.minimum(1,bits[:,:len(R)] @ np.array([c['taticas'] for c in d['controles']])) @ q/q.sum()
    acceptable = feasible & (V >= d['valor_minimo']) & (cov >= d['cobertura_minima'])
    counts = np.random.default_rng(seed).multinomial(shots,prob/prob.sum())
    candidates = np.flatnonzero((counts>0)&feasible)
    best = int(candidates[np.argmax(V[candidates])]) if len(candidates) else None
    gf = greedy(d,R)
    greedy_value = float(gf[:len(R)]@R)
    pfeas = float(prob[feasible].sum())
    return {'Popt':float(prob[optimum].sum()), 'Pground':float(prob[ground].sum()),
            'Pfeas':pfeas,'Paceitavel':float(prob[acceptable].sum()),
            'P90':float(prob[feasible & (V >= .90*opt_value)].sum()),
            'P95':float(prob[feasible & (V >= .95*opt_value)].sum()),
            'Popt_shots':float(counts[optimum].sum()/shots),'Pfeas_shots':float(counts[feasible].sum()/shots),
            'valor_otimo':opt_value,'valor_greedy':greedy_value,
            'valor_condicional_factivel':float(prob[feasible]@V[feasible]/pfeas) if pfeas>0 else None,
            'valor_melhor_medido':float(V[best]) if best is not None else None,
            'gap_melhor_medido':float((opt_value-V[best])/opt_value) if best is not None and opt_value>0 else None,
            'bitstring_melhor_medido':format(best,f'0{len(a)}b') if best is not None else None,
            'otimo_medido':bool(counts[optimum].sum()),
            'controles_melhor_medido':[d['controles'][i]['nome'] for i in range(len(R)) if bits[best,i]] if best is not None else [],
            'aceitabilidade_possivel':bool(acceptable.any())},counts


def run(d,method,p,seed,shots,maxiter,restarts,out,alpha=.2,epsilon=.25,noise=False):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    model,scale,warm = setup(d,method,epsilon)
    R,a,rho,Q,h,J,C,bits,costs,V,E = model
    n=len(a)
    order = np.argsort(E,kind='stable')
    actual_alpha = 1. if method in ('baseline','penalidade_menor') else alpha
    rng=np.random.default_rng(seed)
    histories=[];runs=[]
    started=time.perf_counter()
    for restart in range(restarts):
        initial=np.r_[rng.uniform(0,2*np.pi,p),rng.uniform(0,np.pi,p)]
        trace=[]
        def objective(params):
            prob=abs(estado(params,E,n,scale,warm))**2
            score=cvar(prob,E,actual_alpha,order)/scale
            trace.append(float(score))
            return score
        res=minimize(objective,initial,method='COBYLA',options={'maxiter':maxiter,'rhobeg':1.,'tol':1e-5})
        runs.append({'iniciais':initial.tolist(),'finais':res.x.tolist(),'score':float(res.fun),
                     'avaliacoes':int(res.nfev),'convergiu':bool(res.success),'mensagem':str(res.message)})
        histories.append(trace)
    chosen=min(range(restarts),key=lambda i:runs[i]['score'])
    params=np.array(runs[chosen]['finais'])
    elapsed=time.perf_counter()-started
    prob=abs(estado(params,E,n,scale,warm))**2
    metrics,counts=evaluate(prob,model,d,shots,seed+10000)
    initial_prob=abs(estado(np.zeros(2*p),E,n,scale,warm))**2
    initial_metrics,_=evaluate(initial_prob,model,d,shots,seed+10000)
    qc=circuito(h,J,C,params,E,scale,warm)
    qc_m=qc.copy();qc_m.measure_all()
    tc=transpile(qc_m,basis_gates=['rz','sx','x','cx'],coupling_map=CouplingMap.from_line(n),optimization_level=1,seed_transpiler=seed)
    result={'metodo':method,'controles':len(R),'qubits':n,'p':p,'seed':seed,'shots':shots,
            'orcamento':d['orcamento'],'rho':rho,'escala':scale,'alpha':actual_alpha,
            'epsilon':epsilon if warm is not None else None,'otimizacao_s':elapsed,
            'avaliacoes':sum(r['avaliacoes'] for r in runs),'otimizacoes':runs,
            'score_otimizacao':runs[chosen]['score'],'energia_esperada':float(prob@E),
            'logico':recursos(qc),'transpilado':recursos(tc),
            'Popt_inicial':initial_metrics['Popt'],'Pfeas_inicial':initial_metrics['Pfeas'],
            'warm_probabilidades':warm.tolist() if warm is not None else None,**metrics}
    if noise:
        model_noise=NoiseModel()
        model_noise.add_all_qubit_quantum_error(depolarizing_error(.001,1),['sx','x'])
        model_noise.add_all_qubit_quantum_error(depolarizing_error(.01,2),['cx'])
        model_noise.add_all_qubit_readout_error(ReadoutError([[.98,.02],[.02,.98]]))
        raw=AerSimulator(noise_model=model_noise,method='statevector').run(tc,shots=shots,seed_simulator=seed).result().get_counts()
        nc=np.zeros(len(E),int)
        for b,count in raw.items(): nc[int(b.replace(' ',''),2)]+=count
        # Popt e Pfeas abaixo são frequências sob ruído, não probabilidades exatas.
        nm,_=evaluate(nc/shots,model,d,shots,seed+10000)
        result['ruido']={key:nm[key] for key in ['Popt','Pground','Pfeas','Paceitavel']}
        result['ruido']['shots']=shots
        result['ruido']['modelo']='SX/X:0.001;CX:0.01;leitura:0.02;sem T1/T2'
        np.savez_compressed(out.parent/f'{out.name}_ruido.npz',counts=nc)
    out.mkdir(parents=True,exist_ok=True)
    (out/'resultado.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    (out/'instancia.json').write_text(json.dumps(d,indent=2,ensure_ascii=False))
    np.savez_compressed(out/'distribuicao.npz',probabilidades=prob,counts=counts,energias=E,params=params)
    (out/'circuito_logico.txt').write_text(str(qc_m.draw(output='text')))
    (out/'circuito_transpilado.txt').write_text(str(tc.draw(output='text')))
    with (out/'historico.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['restart','avaliacao','score'])
        for i,trace in enumerate(histories):
            for j,score in enumerate(trace):writer.writerow([i,j+1,score])
    return result


def aggregate(results,out):
    (out/'experimentos.json').write_text(json.dumps(results,indent=2,ensure_ascii=False))
    keys=['metodo','controles','qubits','p','seed','shots','orcamento','rho','escala','alpha','Popt','Pground','Pfeas','Paceitavel','Popt_inicial','valor_otimo','valor_greedy','valor_condicional_factivel','valor_melhor_medido','gap_melhor_medido','otimizacao_s','avaliacoes']
    with (out/'resumo.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows({k:r[k] for k in keys} for r in results)
    methods=[m for m in METHODS if any(r['metodo']==m for r in results)]
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,key in zip(axes,['Popt','Pfeas','otimizacao_s']):
        for m in methods:
            points=[]
            for n in sorted(set(r['controles'] for r in results)):
                vals=[r[key] for r in results if r['metodo']==m and r['controles']==n and r['p']==1]
                if vals:points.append((n,float(np.mean(vals)),float(np.std(vals))))
            if points:
                xs,ys,err=zip(*points);ax.errorbar(xs,ys,yerr=err,marker='o',label=m)
        ax.set(xlabel='Controles',ylabel=key)
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(out/'comparacao.png');plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dados',type=Path,default=Path(__file__).with_name('instancia.json'))
    ap.add_argument('--saida',type=Path,default=Path(__file__).with_name('resultados_melhorias'))
    ap.add_argument('--metodos',nargs='+',choices=METHODS,default=list(METHODS))
    ap.add_argument('--tamanhos',nargs='+',type=int,default=[3,4,5,6,7])
    ap.add_argument('--profundidades',nargs='+',type=int,default=[1])
    ap.add_argument('--seeds',nargs='+',type=int,default=[42,43,44])
    ap.add_argument('--shots',type=int,default=2048)
    ap.add_argument('--maxiter',type=int,default=150)
    ap.add_argument('--restarts',type=int,default=2)
    ap.add_argument('--alpha',type=float,default=.2)
    ap.add_argument('--epsilon',type=float,default=.25)
    ap.add_argument('--ruido',action='store_true',help='Ruído sintético na maior instância para todos os métodos.')
    args=ap.parse_args()
    d=json.loads(args.dados.read_text())
    if not 0<args.alpha<=1 or not 0<args.epsilon<.5:ap.error('alpha em (0,1]; epsilon em (0,0.5).')
    if min(args.shots,args.maxiter,args.restarts,*args.tamanhos,*args.profundidades)<=0:ap.error('Parâmetros positivos necessários.')
    if max(args.tamanhos)>len(d['controles']):ap.error('Tamanho excede a base.')
    args.saida.mkdir(parents=True,exist_ok=True)
    (args.saida/'configuracao.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},indent=2))
    results=[]
    for n in args.tamanhos:
        inst=dict(d);inst['controles']=d['controles'][:n]
        inst['orcamento']=max(1,int(sum(c['custo'] for c in inst['controles'])*d['fracao_orcamento']))
        for p in args.profundidades:
            for seed in args.seeds:
                for method in args.metodos:
                    out=args.saida/f'n{n}_p{p}_seed{seed}_{method}'
                    r=run(inst,method,p,seed,args.shots,args.maxiter,args.restarts,out,args.alpha,args.epsilon,args.ruido and n==max(args.tamanhos))
                    results.append(r)
                    print(f'{n=} {p=} {seed=} {method}: Popt={r["Popt"]:.4f}, Pfeas={r["Pfeas"]:.4f}',flush=True)
    aggregate(results,args.saida)

if __name__=='__main__':main()
