# Redução de decisões: implementação e resultados

Executadas 150 combinações de instância, método e semente. Configuração salva em `configuracao.json`: p=1, 2048 shots, até 100 avaliações COBYLA por início, 2 inicializações por estágio. Referência exata original calculada somente depois da execução de cada método.

| Instância | Método | Qubits pico | Ótima verificada / execuções | Gap médio factível | Execução quântica |
|---|---|---:|---:|---:|---|
| prefixo_3 | qaoa_completo | 5 | 3/3 | 0.00% | sim |
| prefixo_3 | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| prefixo_3 | agrupamento_qaoa | 5 | 3/3 | 0.00% | sim |
| prefixo_3 | rqaoa | 5 | 3/3 | 0.00% | sim |
| prefixo_3 | qiro | 0 | 3/3 | 0.00% | não |
| prefixo_4 | qaoa_completo | 7 | 3/3 | 0.00% | sim |
| prefixo_4 | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| prefixo_4 | agrupamento_qaoa | 6 | 0/3 | 49.48% | sim |
| prefixo_4 | rqaoa | 7 | 1/3 | 33.68% | sim |
| prefixo_4 | qiro | 7 | 3/3 | 0.00% | sim |
| prefixo_5 | qaoa_completo | 8 | 3/3 | 0.00% | sim |
| prefixo_5 | presolve_qaoa | 8 | 3/3 | 0.00% | sim |
| prefixo_5 | agrupamento_qaoa | 5 | 3/3 | 0.00% | sim |
| prefixo_5 | rqaoa | 8 | 1/3 | 17.84% | sim |
| prefixo_5 | qiro | 8 | 1/3 | 1.38% | sim |
| prefixo_6 | qaoa_completo | 9 | 3/3 | 0.00% | sim |
| prefixo_6 | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| prefixo_6 | agrupamento_qaoa | 5 | 0/3 | 31.84% | sim |
| prefixo_6 | rqaoa | 9 | 0/3 | 18.10% | sim |
| prefixo_6 | qiro | 9 | 0/3 | 2.66% | sim |
| prefixo_7 | qaoa_completo | 11 | 3/3 | 0.00% | sim |
| prefixo_7 | presolve_qaoa | 10 | 3/3 | 0.00% | sim |
| prefixo_7 | agrupamento_qaoa | 6 | 0/3 | 1.76% | sim |
| prefixo_7 | rqaoa | 11 | 0/3 | 17.77% | sim |
| prefixo_7 | qiro | 11 | 0/3 | 4.21% | sim |
| custos_mdc5 | qaoa_completo | 13 | 3/3 | 0.00% | sim |
| custos_mdc5 | presolve_qaoa | 11 | 3/3 | 0.00% | sim |
| custos_mdc5 | agrupamento_qaoa | 6 | 0/3 | 1.76% | sim |
| custos_mdc5 | rqaoa | 11 | 0/3 | 17.77% | sim |
| custos_mdc5 | qiro | 11 | 0/3 | 4.21% | sim |
| todos_cabem | qaoa_completo | 12 | 3/3 | 0.00% | sim |
| todos_cabem | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| todos_cabem | agrupamento_qaoa | 6 | 0/3 | 56.20% | sim |
| todos_cabem | rqaoa | 12 | 0/3 | 23.49% | sim |
| todos_cabem | qiro | 12 | 0/3 | 33.62% | sim |
| orcamento_zero | qaoa_completo | 7 | 3/3 | 0.00% | sim |
| orcamento_zero | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| orcamento_zero | agrupamento_qaoa | 3 | 3/3 | 0.00% | sim |
| orcamento_zero | rqaoa | 7 | 1/3 | 0.00% | sim |
| orcamento_zero | qiro | 7 | 3/3 | 0.00% | sim |
| sintetico_8 | qaoa_completo | 11 | 3/3 | 0.00% | sim |
| sintetico_8 | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| sintetico_8 | agrupamento_qaoa | 6 | 0/3 | 26.52% | sim |
| sintetico_8 | rqaoa | 11 | 0/3 | 25.36% | sim |
| sintetico_8 | qiro | 11 | 1/3 | 17.96% | sim |
| sintetico_10 | qaoa_completo | 14 | 1/3 | 0.03% | sim |
| sintetico_10 | presolve_qaoa | 0 | 3/3 | 0.00% | não |
| sintetico_10 | agrupamento_qaoa | 5 | 0/3 | 55.22% | sim |
| sintetico_10 | rqaoa | 14 | 0/3 | 24.86% | sim |
| sintetico_10 | qiro | 14 | 1/3 | 5.84% | sim |

## O que foi implementado

- **Presolve certificado:** portfólios gulosos fornecem limites inferiores; relaxações fracionárias fornecem limites superiores condicionais. Só se fixa uma decisão quando o ramo contrário é inviável ou tem UB estritamente menor que um LB factível. Comparações usam frações exatas dos valores float originais. Empates não são descartados. Cada fixação possui certificado verificável.
- **Núcleo reduzido:** desconta custo e valor dos controles fixados em um, mantém mapa de índices originais, divide custos pelo MDC e substitui orçamento por floor(B/MDC). Limita orçamento à soma dos custos restantes, quando equivalente. Inclui offset constante na energia Ising/QUBO e reconstrói todas as medições para portfólios originais.
- **QAOA no núcleo:** mesma otimização da energia média em todos os métodos QAOA desta comparação, rho=max(R)+0.05 e escala pelo maior coeficiente h/J. Este `qaoa_completo` é uma comparação controlada da redução, não uma reprodução numérica do baseline antigo, que usa outra penalidade/escala. Preservar o ótimo na redução não garante que o circuito o encontre.
- **RQAOA:** correlações Z e ZZ estimadas com shots escolhem fixações zi=±1 ou relações zi=±zj. Substituições conservam exatamente a energia dentro da restrição escolhida e são verificadas por testes exaustivos. A escolha da restrição é heurística e pode descartar o ótimo. O residual Ising é enumerado classicamente até o cutoff. Portfólio bruto inviável e fallback guloso são apresentados separadamente.
- **QIRO inspirado para mochila:** correlações unárias orientam fixações dos controles; uma guarda de orçamento impede escolhas sem completamento factível. Backtracking limitado visita alternativas nos primeiros níveis; DP exata resolve cada núcleo final. Esta é uma adaptação específica, não a reprodução do algoritmo publicado para MIS/MAXSAT. Sem garantia global de otimalidade.
- **Agrupamento por similaridade:** Ward sobre atributos padronizados cria grupos; o melhor valor/custo de cada grupo é o representante e os demais são excluídos. É um experimento heurístico, sem certificado. Custos e valores dos representantes permanecem originais.

## Interpretação e limites

“Ótima verificada” é uma comparação posterior com DP racional no modelo original. A frequência de sucesso entre sementes de RQAOA/QIRO não é Popt de uma distribuição final quântica. Os métodos recursivos fazem múltiplos treinamentos e têm cauda clássica; o total de avaliações e os tempos completos ficam no CSV. Não comparar seus resultados com uma única execução QAOA como se o orçamento computacional fosse idêntico.

Qubits pico é o maior circuito executado. No RQAOA, reduzir o residual não reduz o circuito inicial, e o primeiro estágio ainda precisa caber no simulador. QIRO pode terminar classicamente antes de chamar QAOA se o tamanho já está no cutoff. Presolve que elimina todas as decisões registra zero qubits, sem execução quântica. Isso é contribuição clássica.

As garantias valem para maximizar valor aditivo sob um único orçamento inteiro. Cobertura e valor mínimos continuam pós-avaliados no domínio original; a redução não certifica um ótimo sujeito a essas condições adicionais. Não se usa dominância simples entre controles: dois controles podem ser selecionados juntos.

Ruído sintético foi medido apenas nos métodos QAOA de circuito único em `prefixo_7`, quando `--ruido` foi usado. As correlações recursivas usam shots de circuitos ideais, sem erro de gate/leitura. Não houve reotimização sob ruído. O presolve sem circuito não recebe probabilidades de hardware fictícias.

Instâncias adicionais: custos multiplicados por cinco para testar MDC; orçamento que comporta todos; orçamento zero; controles sintéticos novos com seeds 1008/1010. Dados sintéticos, sem validação organizacional ou evidência de vantagem quântica.

## Reproduzir

```bash
python -m unittest discover -v
python experimentos_reducao.py --ruido
```

## Referências

- [PaPILO: referência conceitual para presolve e reconstrução](https://github.com/scipopt/papilo). O código aqui é próprio e específico de mochila, sem dependência de PaPILO.
- [RQAOA: Bravyi et al.](https://arxiv.org/abs/1910.08980).
- [QIRO: Finžgar et al.](https://arxiv.org/abs/2308.13607).
