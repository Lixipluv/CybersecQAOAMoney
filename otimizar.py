"""Entrada simples para uma instância: presolve certificado é o padrão."""
import argparse
import json
from pathlib import Path
from experimentos_reducao import run,METHODS


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dados',type=Path,default=Path(__file__).with_name('instancia.json'))
    ap.add_argument('--saida',type=Path,default=Path(__file__).with_name('solucao_reduzida'))
    ap.add_argument('--metodo',choices=METHODS,default='presolve_qaoa')
    ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--p',type=int,default=1)
    ap.add_argument('--shots',type=int,default=2048)
    ap.add_argument('--maxiter',type=int,default=100)
    ap.add_argument('--restarts',type=int,default=2)
    ap.add_argument('--cutoff',type=int,default=3)
    ap.add_argument('--grupos',type=int,default=3)
    ap.add_argument('--backtracking',type=int,default=1)
    ap.add_argument('--max-nodes',type=int,default=32)
    ap.add_argument('--ruido',action='store_true')
    args=ap.parse_args()
    if min(args.p,args.shots,args.maxiter,args.restarts,args.cutoff,args.grupos,args.max_nodes)<=0 or args.backtracking<0:
        ap.error('Parâmetros inválidos.')
    args.ruido_todas=True
    result=run(json.loads(args.dados.read_text()),'instancia_usuario',args.metodo,args.seed,args,args.saida)
    print(json.dumps({k:result[k] for k in ['metodo','qubits_modelo_original','qubits_pico','execucao_quantica','reducao_preserva_otimo','otima_verificada','solucao_obtida']},indent=2,ensure_ascii=False))

if __name__=='__main__':main()
