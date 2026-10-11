"""Figuras da evolução experimental, com dados e configurações rastreáveis."""
import csv
import json
from pathlib import Path
from statistics import mean,stdev
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'graficos_evolucao'
COLORS={'inicial':'#65738c','warm':'#9970ab','gm1':'#dd9851','gm3':'#bb6538','completo':'#448bb3','atual':'#1c9275'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white',
                     'savefig.facecolor':'white','axes.labelcolor':'#263247','text.color':'#263247'})


def load(folder):return json.loads((ROOT/folder/'experimentos.json').read_text())

def stats(values):return mean(values),stdev(values) if len(values)>1 else 0.

def finish(fig,name,caption):
    fig.text(.04,.02,caption,fontsize=8,color='#59667b',va='bottom')
    fig.savefig(OUT/f'{name}.png',dpi=200)
    fig.savefig(OUT/f'{name}.svg')
    plt.close(fig)


def collect():
    old=load('resultados');improved=load('resultados_melhorias');gm=load('resultados_mixer');current=load('resultados_reducao')
    groups=[('inicial','Inicial\np=1',[r for r in old if r['controles']==7 and r['p']==1 and r['ambiente']=='ideal']),
            ('warm','Warm\nstart\np=1',[r for r in improved if r['controles']==7 and r['metodo']=='warm_cvar']),
            ('gm1','Mixer\nfactível\np=1',[r for r in gm if r['p']==1]),
            ('gm3','Mixer\nfactível\np=3',[r for r in gm if r['p']==3]),
            ('atual','Atual\npresolve\np=1',[r for r in current if r['instancia']=='prefixo_7' and r['metodo']=='presolve_qaoa'])]
    rows=[]
    for key,label,g in groups:
        if len(g)!=3:raise ValueError(f'Grupo incompleto: {key}')
        for r in g:
            if key=='inicial':
                popt=r['Popt_exata_ideal'];pfeas=r['Pfeas_exata_ideal'];qubits=r['qubits'];cx=r['transpilado']['gates_2q'];depth=r['transpilado']['profundidade']
                value=r['melhor_factivel']['valor'];opt=r['valor_otimo'];evaluations=sum(x['avaliacoes'] for x in r['otimizacoes'])
            elif key=='atual':
                popt=r['probabilidades_ideais']['Popt'];pfeas=r['probabilidades_ideais']['Pfeas'];qubits=r['qubits_pico']
                cx=r['gates_2q_total_circuitos'];depth=r['resultado_metodo']['transpilado']['profundidade']
                value=r['solucao_obtida']['valor'];opt=r['referencia']['valor'];evaluations=r['avaliacoes']
            else:
                popt=r['Popt'];pfeas=r['Pfeas'];qubits=r['qubits'];cx=r['transpilado']['gates_2q'];depth=r['transpilado']['profundidade']
                value=r['valor_melhor_medido'];opt=r['valor_otimo'];evaluations=r['avaliacoes']
            rows.append({'etapa':key,'rotulo':label.replace('\n',' '),'seed':r['seed'],'Popt':popt,'Pfeas':pfeas,
                         'qubits':qubits,'CX':cx,'profundidade':depth,'melhor_valor':value,'valor_otimo':opt,'avaliacoes':evaluations})
    if not np.allclose([r['valor_otimo'] for r in rows],rows[0]['valor_otimo'],rtol=0,atol=1e-10):raise ValueError('Objetivos diferentes.')
    return old,improved,gm,current,groups,rows


def main():
    OUT.mkdir(exist_ok=True)
    old,improved,gm,current,groups,rows=collect()
    with (OUT/'dados_evolucao.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,axes=plt.subplots(1,3,figsize=(13,5.3));xs=np.arange(len(groups));labels=[label for key,label,g in groups]
    for ax,key,title in zip(axes,['Popt','Pfeas','CX'],['Chance do portfólio ótimo','Respeito ao orçamento','Portas de dois qubits']):
        vals=[];errs=[]
        for method,label,g in groups:
            m,s=stats([r[key] for r in rows if r['etapa']==method]);vals.append(m);errs.append(s)
        ax.bar(xs,vals,yerr=errs,capsize=3,color=[COLORS[k] for k,l,g in groups],width=.68)
        ax.set_xticks(xs,labels,fontsize=8);ax.set_title(title,fontsize=11)
        if key!='CX':ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_ylim(0,max(vals)*1.3 if key=='Popt' else 1.17)
        else:ax.set_ylim(0,max(vals)*1.22)
        for x,val,err in zip(xs,vals,errs):ax.text(x,val+err+max(vals)*.035,f'{100*val:.2f}%' if key=='Popt' else f'{100*val:.1f}%' if key=='Pfeas' else f'{val:.0f}',ha='center',fontsize=8)
    fig.suptitle('Da versão inicial ao pipeline atual — 7 controles, orçamento 8',fontsize=15,fontweight='bold')
    fig.subplots_adjust(top=.82,bottom=.22,wspace=.34)
    finish(fig,'01_evolucao_7_controles','Ideal: probabilidades exatas. Médias ± desvio padrão de 3 sementes. Variantes históricas usam até 150 avaliações por início; atual, 100.\nMixer factível usa preparação genérica cara. Comparação descritiva entre estudos, sem orçamento computacional idêntico.')

    # Comparação controlada: mesma política de treinamento e mesmo cenário.
    methods=[('qaoa_completo','Antes do presolve',COLORS['completo']),('presolve_qaoa','Atual com presolve',COLORS['atual'])]
    fig,axes=plt.subplots(1,3,figsize=(12,5.1))
    noise_rows=[]
    for ax,key,title in zip(axes,['Popt','Pfeas','Paceitavel'],['Portfólio ótimo','Dentro do orçamento','Valor e cobertura aceitáveis']):
        for j,(method,label,color) in enumerate(methods):
            g=[r for r in current if r['instancia']=='prefixo_7' and r['metodo']==method]
            vals=[];errs=[]
            for environment,field in [('Ideal','probabilidades_ideais'),('Ruído','frequencias_ruido')]:
                m,s=stats([r[field][key] for r in g]);vals.append(m);errs.append(s)
                noise_rows.append({'metodo':method,'ambiente':environment,'metrica':key,'media':m,'desvio_sementes':s})
            xpos=np.arange(2)+(j-.5)*.32
            ax.bar(xpos,vals,width=.30,yerr=errs,capsize=3,color=color,label=label)
            for x,val,err in zip(xpos,vals,errs):ax.text(x,val+err+(0.001 if key=='Popt' else .02),f'{val*100:.2f}%',ha='center',fontsize=8)
        ax.set_xticks([0,1],['Ideal','Ruído sintético']);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_title(title,fontsize=11)
        ax.set_ylim(0,.025 if key=='Popt' else 1.05)
    axes[0].legend(fontsize=8,loc='upper left');fig.suptitle('Antes × depois do presolve — comparação controlada',fontsize=15,fontweight='bold')
    fig.subplots_adjust(top=.82,bottom=.20,wspace=.34)
    finish(fig,'02_antes_depois_ruido','7 controles, p=1, mesmas 3 sementes, 2 inícios, até 100 avaliações/início e 2048 shots. Ideal exato; ruído: frequências observadas.\nRuído sintético de gates/leitura; sem reotimização sob ruído. Aceitável: valor ≥ 1 e cobertura ponderada ≥ 0,6, além do orçamento.')
    with (OUT/'dados_antes_depois_ruido.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(noise_rows[0]));w.writeheader();w.writerows(noise_rows)

    fig,axes=plt.subplots(1,2,figsize=(11,5.2));sizes=list(range(3,8))
    for label,color,source in [('Inicial',COLORS['inicial'],old),('Atual certificado',COLORS['atual'],current)]:
        q=[];cx=[];err=[]
        for n in sizes:
            if source is old:
                g=[r for r in source if r['controles']==n and r['p']==1 and r['ambiente']=='ideal'];q.append(mean(r['qubits'] for r in g));m,s=stats([r['transpilado']['gates_2q'] for r in g])
            else:
                g=[r for r in source if r['instancia']==f'prefixo_{n}' and r['metodo']=='presolve_qaoa'];q.append(mean(r['qubits_pico'] for r in g));m,s=stats([r['gates_2q_total_circuitos'] for r in g])
            cx.append(m);err.append(s)
        axes[0].plot(sizes,q,'o-',color=color,label=label,linewidth=2)
        axes[1].errorbar(sizes,cx,yerr=err,fmt='o-',color=color,label=label,capsize=3,linewidth=2)
    for ax,title,ylabel in zip(axes,['Qubits necessários','Portas CX transpiladas'],['Qubits pico','CX por circuito']):
        ax.set(title=title,xlabel='Número de controles',ylabel=ylabel,xticks=sizes);ax.legend();ax.grid(axis='y',alpha=.18)
    for n in [3,4,6]:axes[0].annotate('Clássico',xy=(n,0),xytext=(n,.6),ha='center',fontsize=8,color=COLORS['atual'])
    fig.suptitle('Redução de recursos: a decisão pode dispensar o circuito',fontsize=14,fontweight='bold')
    fig.subplots_adjust(top=.81,bottom=.22,wspace=.24)
    finish(fig,'03_recursos_por_tamanho','Zero qubits/CX significa solução integral pelo presolve clássico. Não é sucesso de um circuito quântico.\nEm 5 controles, nenhuma decisão foi eliminada; em 7, uma decisão foi fixada. O número inclui qubits auxiliares de folga.')

    # Mesmo domínio e objetivo para todos os métodos atuais; mostra falhas.
    names=['qaoa_completo','presolve_qaoa','agrupamento_qaoa','rqaoa','qiro']
    labels=['QAOA completo','Presolve certificado','Agrupamento','RQAOA','QIRO adaptado']
    scenarios=list(dict.fromkeys(r['instancia'] for r in current))
    matrix=np.array([[mean(r['otima_verificada'] for r in current if r['instancia']==case and r['metodo']==method) for case in scenarios] for method in names])
    fig,ax=plt.subplots(figsize=(13,5));im=ax.imshow(matrix,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
    ax.set_yticks(range(len(labels)),labels);ax.set_xticks(range(len(scenarios)),[s.replace('prefixo_','').replace('sintetico_','Sintético ').replace('custos_mdc5','MDC ×5').replace('todos_cabem','Todos cabem').replace('orcamento_zero','Orçamento 0') for s in scenarios],rotation=30,ha='right')
    for i in range(len(names)):
        for j,case in enumerate(scenarios):
            g=[r for r in current if r['instancia']==case and r['metodo']==names[i]]
            ax.text(j,i,f'{sum(r["otima_verificada"] for r in g)}/{len(g)}',ha='center',va='center',color='white' if matrix[i,j]>.6 else '#263247')
    fig.colorbar(im,ax=ax,label='Fração de execuções com ótimo verificado',format=PercentFormatter(1))
    ax.set_title('Qual redução conserva bons resultados? — 150 execuções',pad=16,fontsize=14)
    fig.subplots_adjust(left=.17,right=.92,top=.85,bottom=.28)
    finish(fig,'04_sucesso_metodos','Cada célula: execuções ótimas / 3 sementes. A verificação exata ocorre depois do método, no problema original.\nNão é Popt de uma medição. Métodos recursivos têm múltiplos treinamentos e solução clássica do residual; orçamento computacional diferente.')

    # Métrica operacional: melhor valor observado, não energia de QUBOs distintos.
    fig,axes=plt.subplots(1,2,figsize=(12,5.2))
    quality_groups=[('Inicial',COLORS['inicial'],[r for r in rows if r['etapa']=='inicial']),('Warm start',COLORS['warm'],[r for r in rows if r['etapa']=='warm']),('Mixer p=1',COLORS['gm1'],[r for r in rows if r['etapa']=='gm1']),('Atual certificado',COLORS['atual'],[r for r in rows if r['etapa']=='atual'])]
    vals=[mean(r['melhor_valor']/r['valor_otimo'] for r in g) for label,color,g in quality_groups]
    axes[0].barh([label for label,color,g in quality_groups],vals,color=[color for label,color,g in quality_groups]);axes[0].set_xlim(0,1.12);axes[0].xaxis.set_major_formatter(PercentFormatter(1));axes[0].set_title('Melhor valor medido / ótimo clássico')
    for i,value in enumerate(vals):axes[0].text(value+.02,i,f'{100*value:.1f}%',va='center')
    optimal_counts=[];classical_counts=[]
    for method in names:
        g=[r for r in current if r['metodo']==method];optimal_counts.append(sum(r['otima_verificada'] and r['execucao_quantica'] for r in g));classical_counts.append(sum(r['otima_verificada'] and not r['execucao_quantica'] for r in g))
    xs=np.arange(len(names));axes[1].bar(xs,optimal_counts,color=COLORS['completo'],label='Pipeline com circuito');axes[1].bar(xs,classical_counts,bottom=optimal_counts,color=COLORS['atual'],label='Resolvido sem circuito')
    axes[1].set_xticks(xs,['Completo','Presolve','Grupos','RQAOA','QIRO'],rotation=25);axes[1].set_ylim(0,35);axes[1].set_title('Execuções ótimas em 30 casos por método');axes[1].legend(fontsize=8)
    for x,a,b in zip(xs,optimal_counts,classical_counts):axes[1].text(x,a+b+.7,f'{a+b}/30',ha='center',fontsize=9)
    fig.suptitle('Ótimo encontrado e contribuição clássica são medidas diferentes',fontsize=14,fontweight='bold')
    fig.subplots_adjust(left=.14,top=.80,bottom=.23,wspace=.35)
    finish(fig,'05_valor_e_contribuicao_classica','Esquerda: 7 controles, ideal, melhor amostra factível de 2048 shots; todas essas estratégias já atingem o ótimo.\nDireita: 10 cenários × 3 sementes. “Com circuito” inclui componentes clássicos; não atribui o resultado exclusivamente ao QAOA.')

    summary=[]
    for key,label,g in groups:
        subset=[r for r in rows if r['etapa']==key]
        summary.append({'estrategia':label.replace('\n',' '),**{k:stats([r[k] for r in subset])[0] for k in ['Popt','Pfeas','qubits','CX','profundidade','melhor_valor','valor_otimo','avaliacoes']}})
    (OUT/'resumo.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    lines=['# Gráficos da evolução do projeto','',
           'Figuras geradas a partir dos experimentos salvos, sem treinar novamente ou modificar as medições. Há versões PNG em 200 dpi e SVG para ampliação e uso em apresentações. As médias usam três sementes; barras de erro mostram desvio padrão entre sementes, não intervalo de confiança.', '',
           '## Evolução experimental: sete controles','',
           '| Estratégia | Popt ideal | Pfeas ideal | Qubits | CX média |',
           '|---|---:|---:|---:|---:|']
    for r in summary:lines.append(f'| {r["estrategia"]} | {r["Popt"]*100:.2f}% | {r["Pfeas"]*100:.2f}% | {r["qubits"]:.0f} | {r["CX"]:.0f} |')
    lines+=['', '![Evolução experimental](01_evolucao_7_controles.png)', '',
            'A comparação histórica é descritiva: o estudo inicial e o das melhorias usam até 150 avaliações por início; o atual usa até 100. Penalidades, escalas, ansatz e profundidades também diferem. O mixer factível apresenta o melhor comportamento ideal entre estas variantes, mas exige mais portas; o presolve é o padrão atual por preservar o ótimo na redução e economizar recursos.', '',
            '## Antes e depois do presolve: comparação controlada','',
            '![Ideal e ruído](02_antes_depois_ruido.png)', '',
            'Os dois braços usam a mesma configuração de treino e a mesma instância. Ruído disponível neste estudo para sete controles; a versão inicial histórica só simulou ruído em três controles, por isso não foi inventado um ponto inicial de ruído para sete. As frequências sob ruído não são probabilidades exatas. Aceitabilidade é avaliada no portfólio original.', '',
            '## Recursos por tamanho','', '![Recursos](03_recursos_por_tamanho.png)', '',
            'Os zeros nos tamanhos 3, 4 e 6 representam soluções obtidas classicamente pelo presolve. Não são circuitos com zero custo que produziram vantagem quântica.', '',
            '## Qualidade das reduções','', '![Sucesso](04_sucesso_metodos.png)', '',
            'Esta comparação inclui fracassos do agrupamento, RQAOA e QIRO adaptado. Sucesso entre sementes e Popt de uma amostra são grandezas distintas. Os métodos recursivos executam mais treinamentos; a comparação não demonstra superioridade em tempo sob orçamento igual.', '',
            '## Valor obtido e origem do resultado','', '![Valor e contribuição clássica](05_valor_e_contribuicao_classica.png)', '',
            'A versão inicial já encontrava o melhor portfólio dentro das 2048 medições nas instâncias comparadas. As melhorias alteram a frequência de boas medições e os recursos usados. Entre as 30 execuções do pipeline certificado, 21 dispensaram o circuito e nove usaram QAOA no núcleo.', '',
            '## Dados e reprodução','',
            '- `dados_evolucao.csv`: dados por semente da comparação histórica.',
            '- `dados_antes_depois_ruido.csv`: médias e dispersão da comparação controlada.',
            '- `resumo.json`: valores da tabela principal.',
            '- Fontes: `resultados`, `resultados_melhorias`, `resultados_mixer` e `resultados_reducao`, todos no repositório.', '',
            '```bash','python graficos_evolucao.py','```','',
            'As figuras usam o valor original dos controles e métricas de medição. Energias de QUBOs diferentes não são sobrepostas como se fossem a mesma função objetivo. Os dados são sintéticos e não demonstram vantagem quântica.']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n')
    print(f'Criadas 5 figuras PNG/SVG e dados de apoio em {OUT}')

if __name__=='__main__':main()
