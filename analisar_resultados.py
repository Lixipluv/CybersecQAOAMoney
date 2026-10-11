"""Regera o relatório de melhorias a partir de distribuições exportadas."""
import csv
import json
import math
import time
from pathlib import Path
from statistics import mean, stdev
import numpy as np
from projeto import tabela, plt
from mixer_factivel import prepare
from classico import resolver

ROOT=Path(__file__).resolve().parent


def main():
    records=[]
    for folder in ('resultados_melhorias','resultados_mixer'):
        for file in sorted((ROOT/folder).glob('n*/resultado.json')):
            r=json.loads(file.read_text());d=json.loads(file.with_name('instancia.json').read_text())
            prob=np.load(file.with_name('distribuicao.npz'))['probabilidades']
            if r['metodo']=='grover_factivel':
                R,b,cost,V,feasible,initial=prepare(d)
            else:
                R,a,rho,Q,h,J,C,b,cost,V,E=tabela(d,r['rho']);feasible=cost<=d['orcamento']
            opt=float(V[feasible].max())
            r['P90']=float(prob[feasible&(V>=.90*opt)].sum())
            r['P95']=float(prob[feasible&(V>=.95*opt)].sum())
            r['shots_para_95pct_otimo_ideal']=math.ceil(math.log(.05)/math.log1p(-r['Popt'])) if 0<r['Popt']<1 else (1 if r['Popt']>=1 else None)
            file.write_text(json.dumps(r,indent=2,ensure_ascii=False));records.append(r)
    out=ROOT/'resultados_analise';out.mkdir(exist_ok=True)
    (out/'comparacao.json').write_text(json.dumps(records,indent=2,ensure_ascii=False))
    keys=['metodo','controles','qubits','p','seed','Popt','Pfeas','P90','P95','Popt_inicial','valor_otimo','valor_greedy','valor_melhor_medido','otimizacao_s','avaliacoes','shots_para_95pct_otimo_ideal']
    with (out/'comparacao.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k] for k in keys} for r in records)
    largest=max(r['controles'] for r in records)
    seven=[r for r in records if r['controles']==largest]
    groups=[]
    for method,p in dict.fromkeys((r['metodo'],r['p']) for r in seven):
        group=[r for r in seven if r['metodo']==method and r['p']==p]
        groups.append((method,p,group))
    order_methods=['baseline','penalidade_menor','cvar','cvar_escala','warm_cvar','grover_factivel']
    groups.sort(key=lambda g:(order_methods.index(g[0]),g[1]))
    def pct(v):return f'{100*v:.2f}%'
    lines=['# Análise das melhorias: resultados medidos','',
    'Foram testadas 75 combinações de tamanho/semente/estratégia do QAOA com penalidade e nove combinações de profundidade/semente do mixer factível. Cada execução usa duas inicializações e até 150 avaliações por inicialização, com três sementes (42–44) e 2048 shots. As probabilidades ideais são calculadas exatamente pelo statevector; as probabilidades sob ruído são frequências de 2048 shots. O ótimo exato só é usado na avaliação, nunca para escolher parâmetros.', '',
    '## O que estava acontecendo','',
    'O baseline já encontra o valor ótimo nas amostras de todas as 15 execuções de comparação (cinco tamanhos × três sementes). O problema é a baixa frequência do portfólio ótimo, e não a ausência de um portfólio ótimo entre 2048 medições. Na maior instância, uma Popt de aproximadamente 1% torna muito provável observar o ótimo com 2048 shots, mas exige centenas de shots para uma confiança de 95%.', '',
    'O QUBO original tem rho≈5.162 e energia máxima absoluta≈3485.42, contra valor ótimo≈1.85568. Esses números sugeriam testar o condicionamento, mas a comparação mostra que reduzir a penalidade sozinha não resolve a concentração de probabilidade. A energia média considera também os estados com folga inválida; seu mínimo variacional não mede diretamente a qualidade dos portfólios observados.', '',
    '## Comparação na maior instância','',
    '| Estratégia | p | Qubits | Popt ideal média | Pfeas ideal média | P95 ideal média | Popt ruído média | Pfeas ruído média | CX média |',
    '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for method,p,g in groups:
        lines.append(f'| {method} | {p} | {g[0]["qubits"]} | {pct(mean(r["Popt"] for r in g))} | {pct(mean(r["Pfeas"] for r in g))} | {pct(mean(r["P95"] for r in g))} | {pct(mean(r["ruido"]["Popt"] for r in g))} | {pct(mean(r["ruido"]["Pfeas"] for r in g))} | {mean(r["transpilado"]["gates_2q"] for r in g):.0f} |')
    lines+=['', 'P95 é a probabilidade de uma amostra factível obter pelo menos 95% do valor ótimo. Trata-se de qualidade de otimização, não de confiança estatística. Pfeas significa apenas respeitar o orçamento; cobertura e valor mínimos permanecem em Paceitavel. Os dois gráficos incluem desvio padrão entre três sementes, não intervalos de confiança. A contagem CX vem do circuito inteiro transpilado, incluindo a preparação e o mixer.', '',
    '## O que cada mudança realmente oferece','',
    '1. **Penalidade menor:** rho=max(R)+0.05 conserva o ótimo deste modelo de mochila. Foi a alteração mais barata, mas o aumento de Popt em p=1 foi pequeno. Não a tratar como solução milagrosa.',
    '2. **CVaR:** usa a média da cauda inferior de 20% das energias, em vez da energia média completa. Aqui, adicioná-la à penalidade menor não melhorou Popt de modo consistente. O custo CVaR pode ficar plano quando a massa do estado fundamental alcança alpha; portanto, ele não maximiza diretamente Popt. O objetivo difere da energia média pedida no exercício e deve ser apresentado como experimento complementar.',
    '3. **Escala dos coeficientes:** dividir o Hamiltoniano pelo maior coeficiente h/J altera a faixa efetiva de gamma sob a mesma inicialização e os mesmos passos COBYLA. Não muda a expressividade matemática do ansatz; o resultado piorou sob o orçamento de otimização usado. Não selecionar essa escala como nova configuração padrão.',
    '4. **Warm start heurístico:** prepara um produto de qubits com probabilidade 0.75 no bit da heurística valor/custo e 0.25 no outro, incluindo a folga canônica. O mixer é alinhado a esse estado. É uma adaptação heurística, não a reprodução da relaxação contínua do artigo de warm start. No maior tamanho aumenta a frequência ótima, mas reduz a factibilidade; nos tamanhos 3/4 a heurística já é ótima, então o ganho deve ser atribuído principalmente à informação clássica inicial. Em cinco controles houve piora. Popt_inicial está registrada para não atribuir toda a melhora ao circuito.',
    '5. **Mixer factível de Grover:** prepara a superposição uniforme de todos os portfólios dentro do orçamento e mistura dentro desse suporte. Elimina os qubits de folga e as penalidades do circuito experimental; o custo é apenas -sum(R_i x_i). Melhora Popt ideal e garante factibilidade ideal. É um ansatz alternativo, não o circuito QUBO padrão do exercício. A preparação genérica enumera 2^n estados e sintetiza StatePreparation; nesta implementação não escala e é cara em portas. Sob o ruído testado, a vantagem ideal desaparece.', '',
    '## Penalidade menor: derivação usada neste projeto','',
    'Considere custos inteiros positivos, orçamento inteiro não negativo e R_i≥0. Um portfólio inviável excede o orçamento por um inteiro t≥1. Remova controles até torná-lo factível: são necessárias no máximo t remoções, pois cada custo é pelo menos 1. O valor perdido é no máximo t max(R). A energia penalizada do portfólio inviável é no mínimo -V+rho t²; a do portfólio corrigido, com folga canônica, é no máximo -V+t max(R). Se rho>max(R), então rho t²>t max(R), e um inviável não pode ser ótimo. Uma folga incorreta para um portfólio factível também é dominada pela folga canônica. Esta derivação é específica da mochila com uma única restrição e valores não negativos; não generalizar a outras penalidades.', '',
    '## Se a prioridade é obter o melhor portfólio','']
    d=json.loads((ROOT/'instancia.json').read_text())
    d['orcamento']=max(1,int(sum(c['custo'] for c in d['controles'])*d['fracao_orcamento']))
    start=time.perf_counter();dp=resolver(d);runtime=time.perf_counter()-start
    (out/'solucao_classica.json').write_text(json.dumps({**dp,'tempo_s':runtime},indent=2,ensure_ascii=False))
    lines += [f'A programação dinâmica exata (`classico.py`) obtém valor **{dp["valor"]:.6f}**, custo **{dp["custo"]}** dentro do orçamento **{dp["orcamento"]}**. O algoritmo tem O(n B) em transições; usa tuplas para reconstrução do portfólio e é pseudopolinomial no orçamento discretizado. A heurística valor/custo obtém {seven[0]["valor_greedy"]:.6f}, aproximadamente {100*seven[0]["valor_greedy"]/dp["valor"]:.2f}% do ótimo. Para sete controles, a solução exata clássica é a recomendação prática.', '', 'Portfólio ótimo:', '']
    lines += [f'- {c}' for c in dp['controles']]
    baseline=next(g for m,p,g in groups if m=='baseline')
    grover=next(g for m,p,g in groups if m=='grover_factivel' and p==1)
    base_p=mean(r['Popt'] for r in baseline);grover_p=mean(r['Popt'] for r in grover)
    warm=next(g for m,p,g in groups if m=='warm_cvar')
    lines+=['', '## Recomendação para o trabalho','',
    f'Manter a formulação original como baseline e apresentar o mixer factível como comparação: em p=1, Popt ideal média passa de {pct(base_p)} para {pct(grover_p)} ({grover_p/base_p:.2f}×), mas o custo de portas aumenta e o ruído prejudica o resultado. O warm start em sete controles parte de {pct(mean(r["Popt_inicial"] for r in warm))} de frequência ótima e chega a {pct(mean(r["Popt"] for r in warm))}; essa diferença é menor que a comparação com o baseline uniforme pode sugerir.', '',
    'Priorizar Popt, P95, Pfeas, Paceitavel, melhor valor, avaliações clássicas, shots e portas transpiladas em conjunto. Um critério de utilidade deve ser definido antes de ampliar a busca. Não elevar p apenas por expectativa de melhoria: o mixer factível ganhou pouco de p=2 para p=3, com muito mais portas.', '',
    'Próximos testes que justificam investimento: preparação factível usando aritmética reversível/comparadores em vez de StatePreparation genérica; otimização progressiva por profundidade e maior orçamento de busca; reotimização com ruído e modelo calibrado; avaliação em instâncias novas, com custos, ameaças e criticidade diferentes. Não garantir que tais mudanças melhorarão esta ou outras instâncias. Corrigir amostras via heurística pode ser útil, mas os resultados precisam ser rotulados como pós-processamento clássico e comparados com a mesma heurística aplicada a amostras aleatórias.', '',
    'As 84 execuções são explorações em uma mesma base ilustrativa, com poucas sementes; não são validação fora da amostra. O custo CVaR foi calculado com distribuições exatas, e nenhum experimento reotimizou parâmetros sob ruído. Sem evidência de vantagem quântica.', '',
    '## Referências','',
    '- [Barkoutsos et al.: CVaR em otimização variacional](https://arxiv.org/abs/1907.04769).',
    '- [Egger et al.: warm-starting quantum optimization](https://arxiv.org/abs/2009.10095).',
    '- [Hadfield et al.: mixers para restrições](https://arxiv.org/abs/1709.03489).',
    '- [Bärtschi e Eidenbenz: Grover mixers](https://arxiv.org/abs/2006.00354).', '',
    '## Reproduzir','', '```bash', 'python -m unittest discover -v', 'python melhorias.py --ruido', 'python mixer_factivel.py --ruido', 'python analisar_resultados.py', '```']
    (ROOT/'ANALISE_MELHORIAS.md').write_text('\n'.join(lines)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(14,5))
    labels=[m.replace('_','\n')+f'\np={p}' for m,p,g in groups]
    for ax,key in zip(axes,['Popt','Pfeas','gates_2q']):
        vals=[mean(r['transpilado'][key] if key=='gates_2q' else r[key] for r in g) for m,p,g in groups]
        err=[stdev(r['transpilado'][key] if key=='gates_2q' else r[key] for r in g) for m,p,g in groups]
        xs=np.arange(len(groups));ax.bar(xs,vals,yerr=err,color='steelblue',alpha=.85,label='Ideal' if key!='gates_2q' else 'Transpilado')
        if key!='gates_2q':
            noisy=[mean(r['ruido'][key] for r in g) for m,p,g in groups]
            ax.plot(xs,noisy,'o',color='darkorange',label='Ruído (frequência)');ax.set_ylim(0,max(max(vals),max(noisy))*1.25 if key=='Popt' else 1.05)
        ax.set_xticks(xs,labels,fontsize=6);ax.set(ylabel=key,title=f'{largest} controles');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'comparacao_final.png');plt.close(fig)
    print('Relatório gerado:',ROOT/'ANALISE_MELHORIAS.md')

if __name__=='__main__':main()
