"""RQAOA genérico e redução QIRO inspirada, específica de mochila.

Ambas as reduções por correlação são heurísticas. Ótimo global não é consultado.
QIRO aqui usa fixações unárias com guarda de orçamento e backtracking limitado;
não é uma reprodução integral do artigo de QIRO para MIS/MAXSAT.
"""
from fractions import Fraction
import numpy as np
from reducao import core,core_qubo,optimize_ising,exact,incumbent,fractional_bound,reconstruct,rational_json,stage_resources


def eliminate(h,J,C,index,sign,other=None):
    """Substituição exata zi=sign ou zi=sign*zother (Ising triangular)."""
    if sign not in (-1,1) or not 0<=index<len(h) or (other is not None and (other==index or not 0<=other<len(h))):
        raise ValueError('Substituição inválida.')
    h=np.array(h,float,copy=True);A=np.triu(J,1);S=A+A.T;C=float(C)
    keep=[i for i in range(len(h)) if i!=index]
    if other is None:
        C+=sign*h[index]
        for j in keep:h[j]+=sign*S[index,j]
    else:
        h[other]+=sign*h[index];C+=sign*S[index,other]
        for j in keep:
            if j!=other:S[other,j]+=sign*S[index,j];S[j,other]=S[other,j]
    return h[keep],np.triu(S[np.ix_(keep,keep)],1),C


def correlations(prob,bits,shots,seed):
    counts=np.random.default_rng(seed).multinomial(shots,prob/prob.sum())
    z=1-2*bits
    weighted=counts/shots
    single=weighted@z
    pair=z.T@(weighted[:,None]*z)
    return single,pair,counts


def rqaoa(c,v,B,p=1,seed=42,shots=2048,maxiter=150,restarts=2,cutoff=3):
    if cutoff<1:raise ValueError('cutoff positivo necessário.')
    model=core(c,v,B);q=core_qubo(model)
    h,J,C=q['h'],q['J'],q['C'];labels=list(range(len(h)))
    expressions=[(label,1) for label in labels]
    stages=[];evaluations=0;training=0.;rootqubits=len(h)
    while len(h)>cutoff:
        opt=optimize_ising(h,J,C,p,seed+len(stages),maxiter,restarts)
        single,pair,counts=correlations(opt['prob'],opt['bits'],shots,seed+1000+len(stages))
        candidates=[(abs(float(single[i])),i,None,1 if single[i]>=0 else -1) for i in range(len(h))]
        candidates += [(abs(float(pair[i,j])),i,j,1 if pair[i,j]>=0 else -1) for i in range(len(h)) for j in range(i+1,len(h))]
        score,i,j,sign=max(candidates,key=lambda x:x[0])
        eliminated=labels[i];target=labels[j] if j is not None else None
        stages.append({'qubits':len(h),'rotulos':labels.copy(),'variavel_eliminada':eliminated,
                       'relacao':'zi=sign' if j is None else 'zi=sign*zj','alvo':target,'sign':sign,
                       'confianca_correlacao':score,'medias_Z':single.tolist(),
                       'correlacoes_ZZ':pair.tolist(),'shots_correlacao':shots,'otimizacoes':opt['runs'],
                       'avaliacoes':opt['avaliacoes'],'tempo_otimizacao_s':opt['tempo_s'],
                       **stage_resources(h,J,C,opt)})
        expressions=[(target,s*sign) if root==eliminated else (root,s) for root,s in expressions]
        h,J,C=eliminate(h,J,C,i,sign,j);labels.pop(i)
        evaluations+=opt['avaliacoes'];training+=opt['tempo_s']
    # Enumeração apenas do problema residual; sem conhecimento do ótimo original.
    bits=((np.arange(2**len(h))[:,None]>>np.arange(len(h)))&1);z=1-2*bits
    E=z@h+np.einsum('bi,ij,bj->b',z,J,z)+C
    best=int(np.argmin(E));remaining_z={label:int(zi) for label,zi in zip(labels,z[best])}
    full_z=[s if root is None else s*remaining_z[root] for root,s in expressions]
    core_bits=[int((1-zi)//2) for zi in full_z[:len(c)]]
    original=reconstruct(core_bits,model,len(c))
    feasible=sum(ci*bi for ci,bi in zip(c,original))<=B
    fallback=incumbent(c,v,B) if not feasible else None
    return {'bits':original,'factivel':feasible,'fallback_classico_bits':fallback,
            'preserva_otimo':False,'estagios':stages,'qubits_pico':rootqubits,
            'qubits_final':len(h),'avaliacoes':evaluations,'tempo_otimizacao_s':training,
            'estados_residuais_examinados':len(E),'energia_residual':float(E[best]),
            'expressoes_reconstrucao':[{'root':root,'sign':s} for root,s in expressions],
            'estado_Ising_reconstruido':full_z,'rotulos_residuais':labels,'z_residuais':z[best].tolist()}


def qiro(c,v,B,p=1,seed=42,shots=2048,maxiter=150,restarts=2,cutoff=3,backtracking=1,max_nodes=32):
    if cutoff<1 or backtracking<0 or max_nodes<1:raise ValueError('Limites de busca inválidos.')
    stages=[];leaves=[];nodes=0;evaluations=0;training=0.;peak=0
    def visit(fixed,depth):
        nonlocal nodes,evaluations,training,peak
        if nodes>=max_nodes:return
        nodes+=1;node=nodes
        free=[i for i in range(len(c)) if i not in fixed]
        if len(free)<=cutoff:
            bits,value=exact(c,v,B,fixed)
            leaves.append({'bits':bits,'valor':rational_json(value),'fixacoes':{str(i):b for i,b in fixed.items()}})
            return
        model=core(c,v,B,fixed);q=core_qubo(model)
        opt=optimize_ising(q['h'],q['J'],q['C'],p,seed+node,maxiter,restarts)
        singles,pairs,_=correlations(opt['prob'],opt['bits'],shots,seed+1000+node)
        # Somente controles originais são fixados; bits de folga não viram decisões.
        i=max(range(len(free)),key=lambda j:abs(float(singles[j])))
        original=model['free'][i];preferred=int(singles[i]<0)
        branches=[]
        for b in (preferred,1-preferred):
            ub=fractional_bound(c,v,B,{**fixed,original:b})
            if ub is not None:branches.append(b)
        record={'node':node,'profundidade_busca':depth,'qubits':q['n'],'mapa_original':model['free'],
                'fixacoes_anteriores':{str(i):b for i,b in fixed.items()},'variavel':original,
                'bit_preferido':preferred,'confianca':abs(float(singles[i])),
                'ramos_factiveis':branches,'medias_Z':singles[:len(free)].tolist(),
                'correlacoes_ZZ':pairs[:len(free),:len(free)].tolist(),'otimizacoes':opt['runs'],
                'avaliacoes':opt['avaliacoes'],'tempo_otimizacao_s':opt['tempo_s'],'shots_correlacao':shots,
                **stage_resources(q['h'],q['J'],q['C'],opt)}
        stages.append(record);evaluations+=opt['avaliacoes'];training+=opt['tempo_s'];peak=max(peak,q['n'])
        if branches:
            visit({**fixed,original:branches[0]},depth+1)
            if depth<backtracking and len(branches)>1:visit({**fixed,original:branches[1]},depth+1)
    visit({},0)
    # Limite de nós pode impedir chegar a uma folha; fallback explícito.
    if leaves:
        best=max(leaves,key=lambda r:Fraction(r['valor']['numerador'],r['valor']['denominador']))['bits'];fallback=False
    else:best=incumbent(c,v,B);fallback=True
    return {'bits':best,'factivel':True,'preserva_otimo':False,'estagios':stages,'folhas':leaves,
            'qubits_pico':peak,'qubits_final_classico':min(cutoff,len(c)),
            'avaliacoes':evaluations,'tempo_otimizacao_s':training,'nos_explorados':nodes,
            'limite_nos_atingido':nodes>=max_nodes,'backtracking':backtracking,'fallback_classico':fallback}


def cluster_representatives(d,target=3):
    """Agrupamento heurístico; elimina controles sem certificado de dominância."""
    from scipy.cluster.hierarchy import linkage,cut_tree
    from projeto import valores
    R=valores(d);n=len(R)
    if not 1<=target<=n:raise ValueError('Número de grupos inválido.')
    # Sem transformar custo ou valor do objetivo; features apenas definem grupos.
    features=np.array([[c['custo'],c['eficacia']*c['adesao'],c['prioridade'],c['impacto'],
                       *c['taticas'],*c['ativos']] for c in d['controles']],float)
    deviation=features.std(axis=0);deviation[deviation==0]=1
    standardized=(features-features.mean(axis=0))/deviation
    if target==n:labels=np.arange(n)
    else:labels=cut_tree(linkage(standardized,method='ward'),n_clusters=[target]).ravel()
    groups=[[int(i) for i in np.flatnonzero(labels==label)] for label in np.unique(labels)]
    reps=[max(group,key=lambda i:(R[i]/d['controles'][i]['custo'],-i)) for group in groups]
    return {'grupos':groups,'representantes':reps,'fixed':{i:0 for i in range(n) if i not in reps},
            'preserva_otimo':False,'features':'custo,eficacia*adesao,prioridade,impacto,taticas,ativos',
            'metodo':'Ward padronizado + representante valor/custo'}
