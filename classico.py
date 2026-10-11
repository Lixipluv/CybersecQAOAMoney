"""Mochila exata por programação dinâmica O(n B), sem registro de folga."""
import json
from pathlib import Path
from projeto import valores


def resolver(d):
    R=valores(d);B=d['orcamento']
    # Estados de custo exato; tuplas evitam reconstrução com ponteiros mutáveis.
    dp=[None]*(B+1);dp[0]=(0.,())
    for i,c in enumerate(d['controles']):
        cost=c['custo']
        for used in range(B,cost-1,-1):
            previous=dp[used-cost]
            if previous is None:continue
            candidate=(previous[0]+float(R[i]),previous[1]+(i,))
            if dp[used] is None or candidate[0]>dp[used][0]:dp[used]=candidate
    best_cost,best=max(((c,s) for c,s in enumerate(dp) if s is not None),key=lambda p:p[1][0])
    return {'valor':best[0],'custo':best_cost,'orcamento':B,'indices':list(best[1]),
            'controles':[d['controles'][i]['nome'] for i in best[1]],'metodo':'programacao_dinamica_O(nB)'}

if __name__=='__main__':
    print(json.dumps(resolver(json.loads(Path(__file__).with_name('instancia.json').read_text())),indent=2,ensure_ascii=False))
