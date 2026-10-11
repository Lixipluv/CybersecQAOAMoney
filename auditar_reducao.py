"""Verifica certificados, reconstruções RQAOA e integridade dos resultados salvos."""
import argparse
import json
from pathlib import Path
import numpy as np
from reducao import problem,check_certificates,core,core_qubo,exact
from recursivo import eliminate


def audit(directory):
    updated=[];certificates=0;recursions=0
    for path in sorted(directory.glob('*_seed*/resultado.json')):
        record=json.loads(path.read_text());d=json.loads(path.with_name('instancia.json').read_text())
        c,v,B=problem(d);result=record['resultado_metodo'];method=record['metodo']
        if method=='presolve_qaoa':
            red=record['reducao']
            if not check_certificates(c,v,B,red):raise AssertionError(f'Certificado inválido: {path}')
            conditional=exact(c,v,B,{int(i):b for i,b in red['fixed'].items()})[1]
            global_value=exact(c,v,B)[1]
            if conditional!=global_value:raise AssertionError('Presolve alterou o ótimo.')
            certificates+=len(red['certificados'])
        for key in ('distribuicao','distribuicao_ruido'):
            if key in result:
                dist=result[key]
                if not np.isclose(sum(r['prob'] for r in dist),1.,atol=1e-9,rtol=0):raise AssertionError('Probabilidades inválidas.')
                if sum(r['shots'] for r in dist)!=record['shots']:raise AssertionError('Shots inconsistentes.')
                if any(len(r['bits'])!=len(c) or any(b not in (0,1) for b in r['bits']) for r in dist):raise AssertionError('Bits inconsistentes.')
        if method=='rqaoa':
            q=core_qubo(core(c,v,B));h,J,C=q['h'],q['J'],q['C'];labels=list(range(len(h)))
            expressions=[(i,1) for i in labels]
            for stage in result['estagios']:
                i=labels.index(stage['variavel_eliminada']);target=stage['alvo'];sign=stage['sign']
                j=None if target is None else labels.index(target)
                eliminated=labels[i]
                expressions=[(target,s*sign) if root==eliminated else (root,s) for root,s in expressions]
                h,J,C=eliminate(h,J,C,i,sign,j);labels.pop(i)
            bits=((np.arange(2**len(h))[:,None]>>np.arange(len(h)))&1);z=1-2*bits
            E=z@h+np.einsum('bi,ij,bj->b',z,J,z)+C;best=int(np.argmin(E))
            values={label:int(zi) for label,zi in zip(labels,z[best])}
            full=[s if root is None else s*values[root] for root,s in expressions]
            original_energy=np.array(full)@q['h']+np.array(full)@q['J']@np.array(full)+q['C']
            if not np.isclose(original_energy,result['energia_residual'],atol=1e-7,rtol=0):raise AssertionError('Energia RQAOA inconsistente.')
            if [(1-zi)//2 for zi in full[:len(c)]]!=result['bits']:raise AssertionError('Reconstrução inconsistente.')
            result.update({'estado_Ising_reconstruido':full,'rotulos_residuais':labels,'z_residuais':z[best].tolist()})
            recursions+=1
        if 'probabilidades_ideais' in record:
            record['origem_probabilidades']='statevector_ideal' if record['execucao_quantica'] else 'classica_deterministica_sem_circuito'
        record['auditoria_integridade']='verificada'
        path.write_text(json.dumps(record,indent=2,ensure_ascii=False));updated.append(record)
    if not updated:raise ValueError('Nenhuma execução encontrada.')
    # Mantém a mesma ordem do arquivo agregado anterior.
    original=json.loads((directory/'experimentos.json').read_text())
    bykey={(r['instancia'],r['metodo'],r['seed']):r for r in updated}
    ordered=[bykey[(r['instancia'],r['metodo'],r['seed'])] for r in original]
    (directory/'experimentos.json').write_text(json.dumps(ordered,indent=2,ensure_ascii=False))
    summary={'execucoes_auditadas':len(updated),'fixacoes_certificadas_verificadas':certificates,
             'reconstrucoes_RQAOA_verificadas':recursions,'status':'aprovado'}
    (directory/'auditoria.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--saida',type=Path,default=Path(__file__).with_name('resultados_reducao'))
    print(json.dumps(audit(ap.parse_args().saida),indent=2))
