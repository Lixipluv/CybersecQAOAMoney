"""Comparação auditável de redução certificada e heurística de decisões."""
import argparse
import copy
import csv
import json
import time
from pathlib import Path
from statistics import mean
import numpy as np
from projeto import valores,plt
from reducao import problem,presolve,check_certificates,solve_core_qaoa,exact,rational_json
from recursivo import rqaoa,qiro,cluster_representatives

METHODS=['qaoa_completo','presolve_qaoa','agrupamento_qaoa','rqaoa','qiro']


def instances(source,extended=True):
    items=[]
    for n in range(3,min(7,len(source['controles']))+1):
        d=copy.deepcopy(source);d['controles']=d['controles'][:n]
        d['orcamento']=int(sum(c['custo'] for c in d['controles'])*d['fracao_orcamento'])
        items.append((f'prefixo_{n}',d))
    if extended:
        d=copy.deepcopy(source)
        for c in d['controles']:c['custo']*=5
        d['orcamento']=int(sum(c['custo'] for c in d['controles'])*d['fracao_orcamento'])
        items.append(('custos_mdc5',d))
        d=copy.deepcopy(source);d['orcamento']=sum(c['custo'] for c in d['controles']);items.append(('todos_cabem',d))
        d=copy.deepcopy(source);d['orcamento']=0;items.append(('orcamento_zero',d))
        # Duas instâncias novas com parâmetros aleatórios, sementes independentes do treino.
        for n in (8,10):
            rng=np.random.default_rng(1000+n);d=copy.deepcopy(source);d['controles']=[]
            for i in range(n):
                template=copy.deepcopy(source['controles'][i%len(source['controles'])])
                template['nome']=f'Controle sintético {i+1} (n={n})';template['custo']=int(rng.integers(1,6))
                for key in ('eficacia','adesao','prioridade'):template[key]=float(rng.uniform(.4,1.))
                template['impacto']=float(rng.uniform(.05,.6))
                template['taticas']=rng.integers(0,2,len(source['relevancia_taticas'])).tolist()
                template['ativos']=rng.integers(0,2,len(source['ativos'])).tolist()
                d['controles'].append(template)
            d['orcamento']=int(sum(c['custo'] for c in d['controles'])*.4)
            d['descricao']=f'Instância sintética independente, seed={1000+n}; não representa dados organizacionais.'
            items.append((f'sintetico_{n}',d))
    return items


def decode(d,c,v,bits):
    if bits is None:return None
    bits=[int(x) for x in bits]
    if len(bits)!=len(c) or any(x not in (0,1) for x in bits):raise ValueError('Reconstrução inválida.')
    cost=sum(ci*bi for ci,bi in zip(c,bits));value=sum((vi*bi for vi,bi in zip(v,bits)),start=0)
    q=np.array(d['relevancia_taticas'])
    tactics=np.array([ctrl['taticas'] for ctrl in d['controles']])
    coverage=float(np.minimum(1,np.array(bits)@tactics)@q/q.sum())
    return {'bits':bits,'custo':cost,'valor':float(value),'valor_racional':rational_json(value),
            'factivel':cost<=d['orcamento'],'cobertura':coverage,
            'aceitavel':cost<=d['orcamento'] and float(value)>=d['valor_minimo'] and coverage>=d['cobertura_minima'],
            'controles':[ctrl['nome'] for ctrl,b in zip(d['controles'],bits) if b]}


def evaluate_distribution(d,c,v,distribution,optimum):
    popt=pfeas=p95=paceit=0.;shot_opt=shot_feas=total=0
    for record in distribution:
        portfolio=decode(d,c,v,record['bits']);value=sum(vi*bi for vi,bi in zip(v,record['bits']))
        feasible=portfolio['factivel'];optimal=feasible and value==optimum
        popt+=record['prob']*optimal;pfeas+=record['prob']*feasible
        p95+=record['prob']*(feasible and value>=optimum*95/100)
        paceit+=record['prob']*portfolio['aceitavel']
        shot_opt+=record['shots']*optimal;shot_feas+=record['shots']*feasible;total+=record['shots']
    return {'Popt':popt,'Pfeas':pfeas,'P95':p95,'Paceitavel':paceit,
            'Popt_shots':shot_opt/total if total else None,'Pfeas_shots':shot_feas/total if total else None}


def run(d,name,method,seed,args,out):
    c,v,B=problem(d);start=time.perf_counter()
    reduction=None
    kwargs={'p':args.p,'seed':seed,'shots':args.shots,'maxiter':args.maxiter,'restarts':args.restarts}
    if method=='qaoa_completo':
        result=solve_core_qaoa(c,v,B,{},**kwargs,noise=args.ruido and (name=='prefixo_7' or getattr(args,'ruido_todas',False)),simplify=False)
        preserves=True
    elif method=='presolve_qaoa':
        presolve_start=time.perf_counter();reduction=presolve(c,v,B);presolve_time=time.perf_counter()-presolve_start
        if not check_certificates(c,v,B,reduction):raise RuntimeError('Certificados inválidos.')
        result=solve_core_qaoa(c,v,B,reduction['fixed'],**kwargs,noise=args.ruido and (name=='prefixo_7' or getattr(args,'ruido_todas',False)))
        result['tempo_presolve_s']=presolve_time;preserves=True
    elif method=='agrupamento_qaoa':
        reduction=cluster_representatives(d,min(args.grupos,len(c)))
        result=solve_core_qaoa(c,v,B,reduction['fixed'],**kwargs,noise=args.ruido and (name=='prefixo_7' or getattr(args,'ruido_todas',False)));preserves=False
    elif method=='rqaoa':
        result=rqaoa(c,v,B,**kwargs,cutoff=args.cutoff);preserves=False
    else:
        result=qiro(c,v,B,**kwargs,cutoff=args.cutoff,backtracking=args.backtracking,max_nodes=args.max_nodes);preserves=False
    total=time.perf_counter()-start
    # A solução exata original só é calculada depois do método terminar.
    reference_start=time.perf_counter();optimal_bits,optimum=exact(c,v,B);reference_time=time.perf_counter()-reference_start
    decoded=decode(d,c,v,result['bits'])
    obtained=None if decoded is None else sum(vi*bi for vi,bi in zip(v,result['bits']))
    is_optimal=bool(decoded is not None and decoded['factivel'] and obtained==optimum)
    if method in ('rqaoa','qiro'):
        peak=result['qubits_pico'];last=result.get('qubits_final',result.get('qubits_final_classico',0))
        gates=sum(s['transpilado']['gates_2q'] for s in result['estagios'])
        quantum=bool(result['estagios'])
    else:
        peak=result['qubits'];last=result['qubits'];gates=result['transpilado']['gates_2q'];quantum=result['execucao_quantica']
    record={'instancia':name,'metodo':method,'seed':seed,'p':args.p,'shots':args.shots,'controles':len(c),'orcamento':B,
            'qubits_modelo_original':len(c)+B.bit_length(),'qubits_pico':peak,'qubits_final':last,
            'gates_2q_total_circuitos':gates,'avaliacoes':result['avaliacoes'],'tempo_total_metodo_s':total,
            'tempo_otimizacao_s':result['tempo_otimizacao_s'],'tempo_referencia_classica_s':reference_time,
            'execucao_quantica':quantum,'reducao_preserva_otimo':preserves,
            'solucao_obtida':decoded,'otima_verificada':is_optimal,
            'gap_relativo':float((optimum-obtained)/optimum) if decoded and decoded['factivel'] and optimum>0 else (0. if is_optimal else None),
            'referencia':decode(d,c,v,optimal_bits),'reducao':reduction,
            'certificacao_final':'DP racional posterior no modelo original',
            'resultado_metodo':result}
    if 'distribuicao' in result:
        record['origem_probabilidades']='statevector_ideal' if quantum else 'classica_deterministica_sem_circuito'
        record['probabilidades_ideais']=evaluate_distribution(d,c,v,result['distribuicao'],optimum)
    if 'distribuicao_ruido' in result:
        record['frequencias_ruido']=evaluate_distribution(d,c,v,result['distribuicao_ruido'],optimum)
    if result.get('fallback_classico_bits') is not None:
        record['fallback_classico_separado']=decode(d,c,v,result['fallback_classico_bits'])
    out.mkdir(parents=True,exist_ok=True)
    # As fixações passam a chaves string no JSON; check_certificates aceita ambos.
    (out/'resultado.json').write_text(json.dumps(record,indent=2,ensure_ascii=False))
    (out/'instancia.json').write_text(json.dumps(d,indent=2,ensure_ascii=False))
    for key in ('circuito_logico','circuito_transpilado'):
        if key in result:(out/f'{key}.txt').write_text(result[key])
    return record


def report(records,out,args):
    (out/'experimentos.json').write_text(json.dumps(records,indent=2,ensure_ascii=False))
    keys=['instancia','metodo','seed','controles','qubits_modelo_original','qubits_pico','qubits_final','gates_2q_total_circuitos','avaliacoes','tempo_total_metodo_s','execucao_quantica','reducao_preserva_otimo','otima_verificada','gap_relativo']
    with (out/'resumo.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k] for k in keys} for r in records)
    groups=[]
    for name in dict.fromkeys(r['instancia'] for r in records):
        for method in args.metodos:
            g=[r for r in records if r['instancia']==name and r['metodo']==method]
            if g:groups.append((name,method,g))
    lines=['# Redução de decisões: implementação e resultados','',
           f'Executadas {len(records)} combinações de instância, método e semente. Configuração salva em `configuracao.json`: p={args.p}, {args.shots} shots, até {args.maxiter} avaliações COBYLA por início, {args.restarts} inicializações por estágio. Referência exata original calculada somente depois da execução de cada método.', '',
           '| Instância | Método | Qubits pico | Ótima verificada / execuções | Gap médio factível | Execução quântica |',
           '|---|---|---:|---:|---:|---|']
    for name,method,g in groups:
        gaps=[r['gap_relativo'] for r in g if r['gap_relativo'] is not None]
        gap=f'{100*mean(gaps):.2f}%' if gaps else 'sem solução factível'
        quantum='sim' if any(r['execucao_quantica'] for r in g) else 'não'
        lines.append(f'| {name} | {method} | {max(r["qubits_pico"] for r in g)} | {sum(r["otima_verificada"] for r in g)}/{len(g)} | {gap} | {quantum} |')
    lines += ['', '## O que foi implementado','',
              '- **Presolve certificado:** portfólios gulosos fornecem limites inferiores; relaxações fracionárias fornecem limites superiores condicionais. Só se fixa uma decisão quando o ramo contrário é inviável ou tem UB estritamente menor que um LB factível. Comparações usam frações exatas dos valores float originais. Empates não são descartados. Cada fixação possui certificado verificável.',
              '- **Núcleo reduzido:** desconta custo e valor dos controles fixados em um, mantém mapa de índices originais, divide custos pelo MDC e substitui orçamento por floor(B/MDC). Limita orçamento à soma dos custos restantes, quando equivalente. Inclui offset constante na energia Ising/QUBO e reconstrói todas as medições para portfólios originais.',
              '- **QAOA no núcleo:** mesma otimização da energia média em todos os métodos QAOA desta comparação, rho=max(R)+0.05 e escala pelo maior coeficiente h/J. Este `qaoa_completo` é uma comparação controlada da redução, não uma reprodução numérica do baseline antigo, que usa outra penalidade/escala. Preservar o ótimo na redução não garante que o circuito o encontre.',
              '- **RQAOA:** correlações Z e ZZ estimadas com shots escolhem fixações zi=±1 ou relações zi=±zj. Substituições conservam exatamente a energia dentro da restrição escolhida e são verificadas por testes exaustivos. A escolha da restrição é heurística e pode descartar o ótimo. O residual Ising é enumerado classicamente até o cutoff. Portfólio bruto inviável e fallback guloso são apresentados separadamente.',
              '- **QIRO inspirado para mochila:** correlações unárias orientam fixações dos controles; uma guarda de orçamento impede escolhas sem completamento factível. Backtracking limitado visita alternativas nos primeiros níveis; DP exata resolve cada núcleo final. Esta é uma adaptação específica, não a reprodução do algoritmo publicado para MIS/MAXSAT. Sem garantia global de otimalidade.',
              '- **Agrupamento por similaridade:** Ward sobre atributos padronizados cria grupos; o melhor valor/custo de cada grupo é o representante e os demais são excluídos. É um experimento heurístico, sem certificado. Custos e valores dos representantes permanecem originais.', '',
              '## Interpretação e limites','',
              '“Ótima verificada” é uma comparação posterior com DP racional no modelo original. A frequência de sucesso entre sementes de RQAOA/QIRO não é Popt de uma distribuição final quântica. Os métodos recursivos fazem múltiplos treinamentos e têm cauda clássica; o total de avaliações e os tempos completos ficam no CSV. Não comparar seus resultados com uma única execução QAOA como se o orçamento computacional fosse idêntico.', '',
              'Qubits pico é o maior circuito executado. No RQAOA, reduzir o residual não reduz o circuito inicial, e o primeiro estágio ainda precisa caber no simulador. QIRO pode terminar classicamente antes de chamar QAOA se o tamanho já está no cutoff. Presolve que elimina todas as decisões registra zero qubits, sem execução quântica. Isso é contribuição clássica.', '',
              'As garantias valem para maximizar valor aditivo sob um único orçamento inteiro. Cobertura e valor mínimos continuam pós-avaliados no domínio original; a redução não certifica um ótimo sujeito a essas condições adicionais. Não se usa dominância simples entre controles: dois controles podem ser selecionados juntos.', '',
              'Ruído sintético foi medido apenas nos métodos QAOA de circuito único em `prefixo_7`, quando `--ruido` foi usado. As correlações recursivas usam shots de circuitos ideais, sem erro de gate/leitura. Não houve reotimização sob ruído. O presolve sem circuito não recebe probabilidades de hardware fictícias.', '',
              'Instâncias adicionais: custos multiplicados por cinco para testar MDC; orçamento que comporta todos; orçamento zero; controles sintéticos novos com seeds 1008/1010. Dados sintéticos, sem validação organizacional ou evidência de vantagem quântica.', '',
              '## Reproduzir','', '```bash','python -m unittest discover -v','python experimentos_reducao.py --ruido','```','',
              '## Referências','',
              '- [PaPILO: referência conceitual para presolve e reconstrução](https://github.com/scipopt/papilo). O código aqui é próprio e específico de mochila, sem dependência de PaPILO.',
              '- [RQAOA: Bravyi et al.](https://arxiv.org/abs/1910.08980).',
              '- [QIRO: Finžgar et al.](https://arxiv.org/abs/2308.13607).']
    (out/'RELATORIO_REDUCAO.md').write_text('\n'.join(lines)+'\n')
    prefixes=[name for name in dict.fromkeys(r['instancia'] for r in records) if name.startswith('prefixo_')]
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for method in args.metodos:
        x=[];qubits=[];success=[];runtime=[]
        for name in prefixes:
            g=[r for r in records if r['instancia']==name and r['metodo']==method]
            if g:
                x.append(g[0]['controles']);qubits.append(mean(r['qubits_pico'] for r in g))
                success.append(mean(r['otima_verificada'] for r in g));runtime.append(mean(r['tempo_total_metodo_s'] for r in g))
        for ax,y in zip(axes,[qubits,success,runtime]):ax.plot(x,y,marker='o',label=method)
    for ax,label in zip(axes,['Qubits pico','Sucesso ótimo entre sementes','Tempo total por execução (s)']):ax.set(xlabel='Controles',ylabel=label)
    axes[0].legend(fontsize=7);axes[1].set_ylim(-.05,1.05);fig.tight_layout();fig.savefig(out/'comparacao.png');plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dados',type=Path,default=Path(__file__).with_name('instancia.json'))
    ap.add_argument('--saida',type=Path,default=Path(__file__).with_name('resultados_reducao'))
    ap.add_argument('--seeds',type=int,nargs='+',default=[42,43,44])
    ap.add_argument('--metodos',nargs='+',choices=METHODS,default=METHODS)
    ap.add_argument('--p',type=int,default=1)
    ap.add_argument('--shots',type=int,default=2048)
    ap.add_argument('--maxiter',type=int,default=100)
    ap.add_argument('--restarts',type=int,default=2)
    ap.add_argument('--cutoff',type=int,default=3)
    ap.add_argument('--grupos',type=int,default=3)
    ap.add_argument('--backtracking',type=int,default=1)
    ap.add_argument('--max-nodes',type=int,default=32)
    ap.add_argument('--ruido',action='store_true')
    ap.add_argument('--so-prefixos',action='store_true')
    args=ap.parse_args()
    if min(args.p,args.shots,args.maxiter,args.restarts,args.cutoff,args.grupos,args.max_nodes)<=0 or args.backtracking<0:ap.error('Parâmetros inválidos.')
    d=json.loads(args.dados.read_text());items=instances(d,not args.so_prefixos)
    if not items:ap.error('Forneça pelo menos três controles para a varredura.')
    args.saida.mkdir(parents=True,exist_ok=True)
    (args.saida/'configuracao.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},indent=2))
    records=[]
    for name,d in items:
        for seed in args.seeds:
            for method in args.metodos:
                record=run(d,name,method,seed,args,args.saida/f'{name}_seed{seed}_{method}')
                records.append(record)
                print(f'{name} {seed=} {method}: qubits={record["qubits_pico"]}, ótimo={record["otima_verificada"]}',flush=True)
    report(records,args.saida,args)

if __name__=='__main__':main()
