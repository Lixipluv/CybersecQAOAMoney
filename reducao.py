"""Presolve certificado para mochila 0/1, núcleo e QAOA reconstruído.

Frações exatas representam os floats Ri já calculados pelo modelo original.
A certificação é para valor/orçamento; cobertura mínima continua pós-avaliada.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import time
import numpy as np
from scipy.optimize import minimize
from projeto import valores, estado, circuito, recursos


def problem(d):
    return [c['custo'] for c in d['controles']], [Fraction(float(v)) for v in valores(d)], d['orcamento']


def validate(c,v,B):
    if len(c)!=len(v) or not isinstance(B,int) or B<0:
        raise ValueError('Dimensões ou orçamento inválidos.')
    if any(not isinstance(x,int) or x<=0 for x in c) or any(x<0 for x in v):
        raise ValueError('Custos inteiros positivos e valores não negativos necessários.')


def constraints(c,v,B,fixed):
    validate(c,v,B)
    if any(i<0 or i>=len(c) or b not in (0,1) for i,b in fixed.items()):
        raise ValueError('Fixações inválidas.')
    free=[i for i in range(len(c)) if i not in fixed]
    remaining=B-sum(c[i]*b for i,b in fixed.items())
    constant=sum((v[i]*b for i,b in fixed.items()),Fraction(0))
    return free,remaining,constant


def fractional_bound(c,v,B,fixed=None):
    """Limite superior exato da relaxação contínua condicionada às fixações."""
    fixed={} if fixed is None else fixed
    free,remaining,value=constraints(c,v,B,fixed)
    if remaining<0:return None
    for i in sorted(free,key=lambda i:(-v[i]/c[i],i)):
        if remaining>=c[i]:value+=v[i];remaining-=c[i]
        else:
            value+=v[i]*Fraction(remaining,c[i]);break
    return value


def incumbent(c,v,B,fixed=None):
    """Três ordenações gulosas; nenhuma consulta ao ótimo original."""
    fixed={} if fixed is None else fixed
    free,remaining,_=constraints(c,v,B,fixed)
    if remaining<0:raise ValueError('Fixações inviáveis.')
    candidates=[]
    for order in (sorted(free,key=lambda i:(-v[i]/c[i],i)),sorted(free,key=lambda i:(-v[i],i)),sorted(free,key=lambda i:(c[i],i))):
        b=[fixed.get(i,0) for i in range(len(c))];capacity=remaining
        for i in order:
            if c[i]<=capacity:b[i]=1;capacity-=c[i]
        candidates.append(b)
    return max(candidates,key=lambda b:sum((vi*bi for vi,bi in zip(v,b)),Fraction(0)))


def exact(c,v,B,fixed=None):
    """DP racional: referência/certificação posterior ou solução do núcleo final."""
    fixed={} if fixed is None else fixed
    free,remaining,constant=constraints(c,v,B,fixed)
    if remaining<0:raise ValueError('Fixações inviáveis.')
    remaining=min(remaining,sum(c[i] for i in free))
    dp=[None]*(remaining+1);dp[0]=(Fraction(0),())
    for i in free:
        for cost in range(remaining,c[i]-1,-1):
            prev=dp[cost-c[i]]
            if prev is None:continue
            cand=(prev[0]+v[i],prev[1]+(i,))
            if dp[cost] is None or cand[0]>dp[cost][0]:dp[cost]=cand
    best=max((x for x in dp if x is not None),key=lambda x:x[0])
    bits=[fixed.get(i,0) for i in range(len(c))]
    for i in best[1]:bits[i]=1
    return bits,constant+best[0]


def rational_json(value):
    if value is None:return None
    return {'numerador':value.numerator,'denominador':value.denominator,'decimal':float(value)}


def presolve(c,v,B):
    """Fixa apenas se o ramo oposto é inviável ou UB < incumbente factível."""
    fixed={};certificates=[];bound_calls=0
    initial=incumbent(c,v,B)
    while True:
        candidate=incumbent(c,v,B,fixed)
        lower=sum((vi*bi for vi,bi in zip(v,candidate)),Fraction(0))
        changed=False
        for i in range(len(c)):
            if i in fixed:continue
            for excluded in (1,0):
                bound_calls+=1
                ub=fractional_bound(c,v,B,{**fixed,i:excluded})
                if ub is None or ub<lower:
                    chosen=1-excluded
                    certificates.append({'variavel':i,'fixada_em':chosen,'ramo_excluido':excluded,
                                         'fixacoes_anteriores':{str(k):b for k,b in fixed.items()},
                                         'incumbente_bits':candidate,'limite_inferior':rational_json(lower),
                                         'limite_superior_ramo':rational_json(ub),
                                         'motivo':'ramo_inviavel' if ub is None else 'UB_estritamente_menor_que_LB'})
                    fixed[i]=chosen;changed=True;break
            if changed:break
        if not changed:break
    final_incumbent=incumbent(c,v,B,fixed)
    return {'fixed':fixed,'free':[i for i in range(len(c)) if i not in fixed],
            'certificados':certificates,'chamadas_limite':bound_calls,
            'incumbente_inicial':initial,'incumbente_final':final_incumbent,
            'preserva_todos_otimos':True}


def check_certificates(c,v,B,result):
    """Reexecuta e verifica as desigualdades exatas e cada incumbente."""
    fixed={}
    for cert in result['certificados']:
        if {int(k):b for k,b in cert['fixacoes_anteriores'].items()}!=fixed:return False
        i=cert['variavel'];excluded=cert['ramo_excluido'];chosen=cert['fixada_em']
        if i in fixed or excluded not in (0,1) or chosen!=1-excluded:return False
        b=cert['incumbente_bits']
        if len(b)!=len(c) or any(x not in (0,1) for x in b):return False
        if any(b[k]!=x for k,x in fixed.items()) or sum(ci*bi for ci,bi in zip(c,b))>B:return False
        lower=sum((vi*bi for vi,bi in zip(v,b)),Fraction(0))
        ub=fractional_bound(c,v,B,{**fixed,i:excluded})
        if rational_json(lower)!=cert['limite_inferior'] or rational_json(ub)!=cert['limite_superior_ramo']:return False
        if ub is not None and ub>=lower:return False
        fixed[i]=chosen
    return fixed=={int(k):b for k,b in result['fixed'].items()} and result['free']==[i for i in range(len(c)) if i not in fixed]


def core(c,v,B,fixed=None,simplify=True):
    """Orçamento residual, offset e redução exata por MDC dos custos."""
    fixed={} if fixed is None else dict(fixed)
    free,remaining,constant=constraints(c,v,B,fixed)
    if remaining<0:raise ValueError('Núcleo inviável.')
    divisor=(reduce(gcd,(c[i] for i in free),0) or 1) if simplify else 1
    scaled_costs=[c[i]//divisor for i in free]
    capacity=min(remaining//divisor,sum(scaled_costs)) if simplify else remaining
    return {'free':free,'fixed':fixed,'costs':scaled_costs,'values':[v[i] for i in free],
            'budget':capacity,'offset_value':constant,'gcd':divisor,'remaining_budget_original':remaining}


def reconstruct(reduced_bits,model,n):
    if len(reduced_bits)!=len(model['free']):raise ValueError('Dimensão do núcleo incorreta.')
    if any(b not in (0,1) for b in reduced_bits):raise ValueError('Bits inválidos.')
    bits=[model['fixed'].get(i,0) for i in range(n)]
    for i,b in zip(model['free'],reduced_bits):bits[i]=int(b)
    return bits


def core_qubo(model):
    """QUBO do núcleo com energia original reconstruída (offset constante)."""
    k=len(model['free']);B=model['budget'];m=B.bit_length()
    if k==0:raise ValueError('Núcleo vazio: nenhuma execução quântica necessária.')
    R=np.array([float(x) for x in model['values']])
    a=np.array(model['costs']+[2**j for j in range(m)],float)
    rho=float(R.max()+.05)
    linear=rho*(a*a-2*B*a);linear[:k]-=R
    Q=np.diag(linear)
    for i in range(len(a)):
        for j in range(i+1,len(a)):Q[i,j]=2*rho*a[i]*a[j]
    J=np.triu(Q,1)/4;h=-linear/2-J.sum(axis=0)-J.sum(axis=1)
    C=rho*B*B+linear.sum()/2+J.sum()-float(model['offset_value'])
    n=len(a)
    if n>20:raise ValueError('Simulação limitada a 20 qubits.')
    b=((np.arange(2**n)[:,None]>>np.arange(n))&1).astype(int)
    E=-(b[:,:k]@R)-float(model['offset_value'])+rho*(b@a-B)**2
    return {'Q':Q,'h':h,'J':J,'C':C,'rho':rho,'bits':b,'E':E,'n':n,'k':k,
            'escala':max(float(np.max(abs(h))),float(np.max(abs(J))),1.)}


def optimize_ising(h,J,C,p,seed,maxiter,restarts,scale=None):
    """Treino ideal sem usar ótimo/bitstring de referência como objetivo."""
    n=len(h)
    if not 1<=n<=20:raise ValueError('Exige 1 a 20 qubits.')
    if p<1 or maxiter<1 or restarts<1:raise ValueError('Parâmetros positivos necessários.')
    bits=((np.arange(2**n)[:,None]>>np.arange(n))&1)
    z=1-2*bits
    E=z@h+np.einsum('bi,ij,bj->b',z,J,z)+C
    scale=max(float(np.max(abs(h))),float(np.max(abs(J))),1.) if scale is None else scale
    rng=np.random.default_rng(seed);runs=[];start=time.perf_counter()
    for restart in range(restarts):
        initial=np.r_[rng.uniform(0,2*np.pi,p),rng.uniform(0,np.pi,p)]
        def objective(x):return float((abs(estado(x,E,n,scale))**2)@E)/scale
        res=minimize(objective,initial,method='COBYLA',options={'maxiter':maxiter,'rhobeg':1.,'tol':1e-5})
        runs.append({'iniciais':initial.tolist(),'finais':res.x.tolist(),'energia_final':float(res.fun*scale),
                     'avaliacoes':int(res.nfev),'convergiu':bool(res.success),'mensagem':str(res.message)})
    best=min(runs,key=lambda r:r['energia_final']);params=np.array(best['finais'])
    prob=abs(estado(params,E,n,scale))**2
    return {'prob':prob,'params':params,'E':E,'bits':bits,'runs':runs,'scale':scale,
            'avaliacoes':sum(r['avaliacoes'] for r in runs),'tempo_s':time.perf_counter()-start}


def solve_core_qaoa(c,v,B,fixed,p=1,seed=42,shots=2048,maxiter=150,restarts=2,noise=False,simplify=True):
    """Reconstrói todas as medições no domínio original; sem reparo oculto."""
    model=core(c,v,B,fixed,simplify);k=len(model['free'])
    result={'free':model['free'],'fixed':model['fixed'],'gcd_custos':model['gcd'],
            'orcamento_nucleo':model['budget'],'offset_valor':float(model['offset_value'])}
    if not k:
        bits=reconstruct([],model,len(c))
        return {**result,'bits':bits,'qubits':0,'avaliacoes':0,'tempo_otimizacao_s':0.,
                'execucao_quantica':False,'distribuicao':[{'bits':bits,'prob':1.,'shots':shots}],
                'logico':{'qubits':0,'profundidade':0,'gates':0,'gates_2q':0},
                'transpilado':{'qubits':0,'profundidade':0,'gates':0,'gates_2q':0}}
    q=core_qubo(model)
    opt=optimize_ising(q['h'],q['J'],q['C'],p,seed,maxiter,restarts)
    qc=circuito(q['h'],q['J'],q['C'],opt['params'],opt['E'],opt['scale'])
    measured=qc.copy();measured.measure_all()
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    tc=transpile(measured,basis_gates=['rz','sx','x','cx'],coupling_map=CouplingMap.from_line(q['n']),optimization_level=1,seed_transpiler=seed)
    prob=opt['prob'];rng=np.random.default_rng(seed+10000);counts=rng.multinomial(shots,prob/prob.sum())
    def decode_distribution(probabilities,frequencies):
        merged={}
        for i,b in enumerate(q['bits']):
            full=tuple(reconstruct(b[:k],model,len(c)))
            rec=merged.setdefault(full,{'bits':list(full),'prob':0.,'shots':0})
            rec['prob']+=float(probabilities[i]);rec['shots']+=int(frequencies[i])
        return list(merged.values())
    distribution=decode_distribution(prob,counts)
    measured_candidates=[r['bits'] for r in distribution if r['shots'] and sum(ci*bi for ci,bi in zip(c,r['bits']))<=B]
    best=max(measured_candidates,key=lambda b:sum(vi*bi for vi,bi in zip(v,b))) if measured_candidates else None
    result.update({'bits':best,'qubits':q['n'],'avaliacoes':opt['avaliacoes'],'tempo_otimizacao_s':opt['tempo_s'],
                   'execucao_quantica':True,'distribuicao':distribution,'logico':recursos(qc),
                   'transpilado':recursos(tc),'otimizacoes':opt['runs'],'rho':q['rho'],'escala':opt['scale'],
                   'circuito_logico':str(measured.draw(output='text')),
                   'circuito_transpilado':str(tc.draw(output='text'))})
    if noise:
        from qiskit_aer import AerSimulator
        from qiskit_aer.noise import NoiseModel,depolarizing_error,ReadoutError
        nm=NoiseModel();nm.add_all_qubit_quantum_error(depolarizing_error(.001,1),['sx','x'])
        nm.add_all_qubit_quantum_error(depolarizing_error(.01,2),['cx'])
        nm.add_all_qubit_readout_error(ReadoutError([[.98,.02],[.02,.98]]))
        raw=AerSimulator(noise_model=nm,method='statevector').run(tc,shots=shots,seed_simulator=seed).result().get_counts()
        noisy=np.zeros(len(prob),int)
        for b,nshots in raw.items():noisy[int(b.replace(' ',''),2)]+=nshots
        result['distribuicao_ruido']=decode_distribution(noisy/shots,noisy)
    return result


def stage_resources(h,J,C,opt):
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap
    qc=circuito(h,J,C,opt['params'],opt['E'],opt['scale'])
    tc=transpile(qc,basis_gates=['rz','sx','x','cx'],coupling_map=CouplingMap.from_line(len(h)),optimization_level=1,seed_transpiler=42)
    return {'logico':recursos(qc),'transpilado':recursos(tc)}
