"""GM-QAOA exploratório: suporte factível e preparação genérica exponencial.

A enumeração é usada para preparar a superposição factível e para avaliar.
Não há consulta ao ótimo no otimizador. Não é uma implementação escalável.
"""
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from projeto import valores, recursos, plt
from melhorias import cvar, greedy


def prepare(d):
    R=valores(d);n=len(R)
    if n>14:raise ValueError('Preparação genérica limitada a 14 controles.')
    b=((np.arange(2**n)[:,None]>>np.arange(n))&1)
    cost=b@np.array([c['custo'] for c in d['controles']]);V=b@R
    feasible=cost<=d['orcamento']
    # Orçamento, não ótimo/valor, determina o suporte inicial.
    initial=feasible.astype(complex)/np.sqrt(feasible.sum())
    return R,b,cost,V,feasible,initial


def state(params,E,initial,scale):
    psi=initial.copy();p=len(params)//2
    for gamma,beta in zip(params[:p],params[p:]):
        psi*=np.exp(-1j*gamma*E/scale)
        psi+=(np.exp(-1j*beta)-1)*np.vdot(initial,psi)*initial
    return psi


def build(R,initial,params,scale):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import StatePreparation
    n=len(R);p=len(params)//2
    A=StatePreparation(initial)
    qc=QuantumCircuit(n);qc.append(A,range(n))
    for gamma,beta in zip(params[:p],params[p:]):
        # -R_i x_i = -R_i/2 + R_i Z_i/2.
        for i in range(n):qc.rz(gamma*R[i]/scale,i)
        qc.global_phase+=gamma*R.sum()/(2*scale)
        # exp(-i beta |F><F|) = A diag(exp(-i beta),1,...) A†.
        qc.append(A.inverse(),range(n));qc.x(range(n))
        if n==1:qc.p(-beta,0)
        else:qc.mcp(-beta,list(range(n-1)),n-1)
        qc.x(range(n));qc.append(A,range(n))
    return qc


def run(d,p,seed,shots,maxiter,restarts,out,alpha=.2,noise=False):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    preparation_start=time.perf_counter()
    R,b,cost,V,feasible,initial=prepare(d)
    prep_time=time.perf_counter()-preparation_start
    E=-V;scale=max(float(R.sum()),1.);order=np.argsort(E,kind='stable')
    rng=np.random.default_rng(seed);runs=[];histories=[]
    started=time.perf_counter()
    for rep in range(restarts):
        x0=np.r_[rng.uniform(0,2*np.pi,p),rng.uniform(0,2*np.pi,p)]
        trace=[]
        def objective(params):
            prob=abs(state(params,E,initial,scale))**2
            score=cvar(prob,E,alpha,order)/scale
            trace.append(score);return score
        res=minimize(objective,x0,method='COBYLA',options={'maxiter':maxiter,'rhobeg':1.,'tol':1e-5})
        runs.append({'iniciais':x0.tolist(),'finais':res.x.tolist(),'score':float(res.fun),
                     'avaliacoes':int(res.nfev),'convergiu':bool(res.success),'mensagem':str(res.message)})
        histories.append(trace)
    chosen=min(range(restarts),key=lambda i:runs[i]['score']);params=np.array(runs[chosen]['finais'])
    elapsed=time.perf_counter()-started
    prob=abs(state(params,E,initial,scale))**2
    optimum=feasible&np.isclose(V,V[feasible].max(),atol=1e-10,rtol=0)
    q=np.array(d['relevancia_taticas'])
    cov=np.minimum(1,b@np.array([c['taticas'] for c in d['controles']]))@q/q.sum()
    acceptable=feasible&(V>=d['valor_minimo'])&(cov>=d['cobertura_minima'])
    counts=rng.multinomial(shots,prob/prob.sum());cand=np.flatnonzero((counts>0)&feasible)
    best=int(cand[np.argmax(V[cand])])
    def metrics(distribution):
        return {'Popt':float(distribution[optimum].sum()),'Pfeas':float(distribution[feasible].sum()),
                'Paceitavel':float(distribution[acceptable].sum()),'valor_esperado':float(distribution@V),
                'P90':float(distribution[feasible & (V >= .90*V[feasible].max())].sum()),
                'P95':float(distribution[feasible & (V >= .95*V[feasible].max())].sum())}
    qc=build(R,initial,params,scale);qm=qc.copy();qm.measure_all()
    transpile_start=time.perf_counter()
    tc=transpile(qm,basis_gates=['rz','sx','x','cx'],coupling_map=CouplingMap.from_line(len(R)),optimization_level=1,seed_transpiler=seed)
    transpile_time=time.perf_counter()-transpile_start
    result={'metodo':'grover_factivel','controles':len(R),'qubits':len(R),'p':p,'seed':seed,'shots':shots,
            'orcamento':d['orcamento'],'alpha':alpha,'escala':scale,'otimizacao_s':elapsed,
            'preparacao_classica_s':prep_time,'transpilacao_s':transpile_time,
            'avaliacoes':sum(r['avaliacoes'] for r in runs),'otimizacoes':runs,
            'estados_factiveis':int(feasible.sum()),'Popt_inicial':float(abs(initial[optimum])@abs(initial[optimum])),
            'logico':recursos(qc),'transpilado':recursos(tc),
            'valor_otimo':float(V[feasible].max()),'valor_greedy':float(greedy(d,R)[:len(R)]@R),
            'valor_melhor_medido':float(V[best]),'bitstring_melhor_medido':format(best,f'0{len(R)}b'),
            'Popt_shots':float(counts[optimum].sum()/shots),'Pfeas_shots':float(counts[feasible].sum()/shots),
            'aceitabilidade_possivel':bool(acceptable.any()),
            'preparacao':'StatePreparation genérica a partir de enumeração 2^n',**metrics(prob)}
    if noise:
        noise_model=NoiseModel()
        noise_model.add_all_qubit_quantum_error(depolarizing_error(.001,1),['sx','x'])
        noise_model.add_all_qubit_quantum_error(depolarizing_error(.01,2),['cx'])
        noise_model.add_all_qubit_readout_error(ReadoutError([[.98,.02],[.02,.98]]))
        raw=AerSimulator(noise_model=noise_model,method='statevector').run(tc,shots=shots,seed_simulator=seed).result().get_counts()
        noisy=np.zeros(len(E),int)
        for bits,count in raw.items():noisy[int(bits.replace(' ',''),2)]+=count
        result['ruido']={**metrics(noisy/shots),'shots':shots,'modelo':'SX/X:0.001;CX:0.01;leitura:0.02;sem T1/T2'}
    out.mkdir(parents=True,exist_ok=True)
    (out/'resultado.json').write_text(json.dumps(result,indent=2,ensure_ascii=False))
    (out/'instancia.json').write_text(json.dumps(d,indent=2,ensure_ascii=False))
    np.savez_compressed(out/'distribuicao.npz',probabilidades=prob,counts=counts,energias=E,params=params)
    (out/'circuito_logico.txt').write_text(str(qm.draw(output='text')))
    (out/'circuito_transpilado.txt').write_text(str(tc.draw(output='text')))
    with (out/'historico.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['restart','avaliacao','score'])
        for i,trace in enumerate(histories):
            for j,score in enumerate(trace):writer.writerow([i,j+1,score])
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dados',type=Path,default=Path(__file__).with_name('instancia.json'))
    ap.add_argument('--saida',type=Path,default=Path(__file__).with_name('resultados_mixer'))
    ap.add_argument('--tamanho',type=int,default=7)
    ap.add_argument('--profundidades',type=int,nargs='+',default=[1,2,3])
    ap.add_argument('--seeds',type=int,nargs='+',default=[42,43,44])
    ap.add_argument('--shots',type=int,default=2048)
    ap.add_argument('--maxiter',type=int,default=150)
    ap.add_argument('--restarts',type=int,default=2)
    ap.add_argument('--alpha',type=float,default=.2)
    ap.add_argument('--ruido',action='store_true')
    args=ap.parse_args();d=json.loads(args.dados.read_text())
    if not 0<args.alpha<=1 or min(args.tamanho,args.shots,args.maxiter,args.restarts,*args.profundidades)<=0:
        ap.error('Parâmetros inválidos.')
    if args.tamanho>len(d['controles']):ap.error('Tamanho excede a base.')
    d['controles']=d['controles'][:args.tamanho]
    d['orcamento']=max(1,int(sum(c['custo'] for c in d['controles'])*d['fracao_orcamento']))
    args.saida.mkdir(parents=True,exist_ok=True)
    (args.saida/'configuracao.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},indent=2))
    results=[]
    for p in args.profundidades:
        for seed in args.seeds:
            r=run(d,p,seed,args.shots,args.maxiter,args.restarts,args.saida/f'n{args.tamanho}_p{p}_seed{seed}',args.alpha,args.ruido)
            results.append(r);print(f'{p=} {seed=}: Popt={r["Popt"]:.4f}, Pfeas={r["Pfeas"]:.4f}, CX={r["transpilado"]["gates_2q"]}',flush=True)
    (args.saida/'experimentos.json').write_text(json.dumps(results,indent=2,ensure_ascii=False))
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,key in zip(axes,['Popt','Pfeas','gates_2q']):
        for p in args.profundidades:
            vals=[r['transpilado'][key] if key=='gates_2q' else r[key] for r in results if r['p']==p]
            ax.errorbar(p,np.mean(vals),yerr=np.std(vals),fmt='o',color='steelblue')
        ax.set(xlabel='p',ylabel=key)
    fig.tight_layout();fig.savefig(args.saida/'comparacao.png');plt.close(fig)

if __name__=='__main__':main()
