# Análise das melhorias: resultados medidos

Foram testadas 75 combinações de tamanho/semente/estratégia do QAOA com penalidade e nove combinações de profundidade/semente do mixer factível. Cada execução usa duas inicializações e até 150 avaliações por inicialização, com três sementes (42–44) e 2048 shots. As probabilidades ideais são calculadas exatamente pelo statevector; as probabilidades sob ruído são frequências de 2048 shots. O ótimo exato só é usado na avaliação, nunca para escolher parâmetros.

## O que estava acontecendo

O baseline já encontra o valor ótimo nas amostras de todas as 15 execuções de comparação (cinco tamanhos × três sementes). O problema é a baixa frequência do portfólio ótimo, e não a ausência de um portfólio ótimo entre 2048 medições. Na maior instância, uma Popt de aproximadamente 1% torna muito provável observar o ótimo com 2048 shots, mas exige centenas de shots para uma confiança de 95%.

O QUBO original tem rho≈5.162 e energia máxima absoluta≈3485.42, contra valor ótimo≈1.85568. Esses números sugeriam testar o condicionamento, mas a comparação mostra que reduzir a penalidade sozinha não resolve a concentração de probabilidade. A energia média considera também os estados com folga inválida; seu mínimo variacional não mede diretamente a qualidade dos portfólios observados.

## Comparação na maior instância

| Estratégia | p | Qubits | Popt ideal média | Pfeas ideal média | P95 ideal média | Popt ruído média | Pfeas ruído média | CX média |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 1 | 11 | 0.98% | 76.80% | 11.00% | 1.04% | 61.20% | 245 |
| penalidade_menor | 1 | 11 | 0.98% | 77.19% | 11.00% | 1.07% | 61.46% | 245 |
| cvar | 1 | 11 | 0.97% | 76.90% | 10.82% | 1.09% | 61.49% | 245 |
| cvar_escala | 1 | 11 | 0.75% | 56.99% | 9.17% | 0.94% | 49.69% | 245 |
| warm_cvar | 1 | 11 | 1.74% | 68.67% | 13.39% | 1.50% | 54.28% | 245 |
| grover_factivel | 1 | 7 | 4.28% | 100.00% | 50.99% | 0.78% | 42.27% | 1298 |
| grover_factivel | 2 | 7 | 5.61% | 100.00% | 62.19% | 0.67% | 40.76% | 2269 |
| grover_factivel | 3 | 7 | 5.86% | 100.00% | 63.52% | 0.88% | 39.55% | 3244 |

P95 é a probabilidade de uma amostra factível obter pelo menos 95% do valor ótimo. Trata-se de qualidade de otimização, não de confiança estatística. Pfeas significa apenas respeitar o orçamento; cobertura e valor mínimos permanecem em Paceitavel. Os dois gráficos incluem desvio padrão entre três sementes, não intervalos de confiança. A contagem CX vem do circuito inteiro transpilado, incluindo a preparação e o mixer.

## O que cada mudança realmente oferece

1. **Penalidade menor:** rho=max(R)+0.05 conserva o ótimo deste modelo de mochila. Foi a alteração mais barata, mas o aumento de Popt em p=1 foi pequeno. Não a tratar como solução milagrosa.
2. **CVaR:** usa a média da cauda inferior de 20% das energias, em vez da energia média completa. Aqui, adicioná-la à penalidade menor não melhorou Popt de modo consistente. O custo CVaR pode ficar plano quando a massa do estado fundamental alcança alpha; portanto, ele não maximiza diretamente Popt. O objetivo difere da energia média pedida no exercício e deve ser apresentado como experimento complementar.
3. **Escala dos coeficientes:** dividir o Hamiltoniano pelo maior coeficiente h/J altera a faixa efetiva de gamma sob a mesma inicialização e os mesmos passos COBYLA. Não muda a expressividade matemática do ansatz; o resultado piorou sob o orçamento de otimização usado. Não selecionar essa escala como nova configuração padrão.
4. **Warm start heurístico:** prepara um produto de qubits com probabilidade 0.75 no bit da heurística valor/custo e 0.25 no outro, incluindo a folga canônica. O mixer é alinhado a esse estado. É uma adaptação heurística, não a reprodução da relaxação contínua do artigo de warm start. No maior tamanho aumenta a frequência ótima, mas reduz a factibilidade; nos tamanhos 3/4 a heurística já é ótima, então o ganho deve ser atribuído principalmente à informação clássica inicial. Em cinco controles houve piora. Popt_inicial está registrada para não atribuir toda a melhora ao circuito.
5. **Mixer factível de Grover:** prepara a superposição uniforme de todos os portfólios dentro do orçamento e mistura dentro desse suporte. Elimina os qubits de folga e as penalidades do circuito experimental; o custo é apenas -sum(R_i x_i). Melhora Popt ideal e garante factibilidade ideal. É um ansatz alternativo, não o circuito QUBO padrão do exercício. A preparação genérica enumera 2^n estados e sintetiza StatePreparation; nesta implementação não escala e é cara em portas. Sob o ruído testado, a vantagem ideal desaparece.

## Penalidade menor: derivação usada neste projeto

Considere custos inteiros positivos, orçamento inteiro não negativo e R_i≥0. Um portfólio inviável excede o orçamento por um inteiro t≥1. Remova controles até torná-lo factível: são necessárias no máximo t remoções, pois cada custo é pelo menos 1. O valor perdido é no máximo t max(R). A energia penalizada do portfólio inviável é no mínimo -V+rho t²; a do portfólio corrigido, com folga canônica, é no máximo -V+t max(R). Se rho>max(R), então rho t²>t max(R), e um inviável não pode ser ótimo. Uma folga incorreta para um portfólio factível também é dominada pela folga canônica. Esta derivação é específica da mochila com uma única restrição e valores não negativos; não generalizar a outras penalidades.

## Se a prioridade é obter o melhor portfólio

A programação dinâmica exata (`classico.py`) obtém valor **1.855676**, custo **7** dentro do orçamento **8**. O algoritmo tem O(n B) em transições; usa tuplas para reconstrução do portfólio e é pseudopolinomial no orçamento discretizado. A heurística valor/custo obtém 1.761023, aproximadamente 94.90% do ótimo. Para sete controles, a solução exata clássica é a recomendação prática.

Portfólio ótimo:

- CIS-1: Inventário de hardware
- CIS-2: Inventário de software
- CIS-10: Defesa contra malware

## Recomendação para o trabalho

Manter a formulação original como baseline e apresentar o mixer factível como comparação: em p=1, Popt ideal média passa de 0.98% para 4.28% (4.37×), mas o custo de portas aumenta e o ruído prejudica o resultado. O warm start em sete controles parte de 1.48% de frequência ótima e chega a 1.74%; essa diferença é menor que a comparação com o baseline uniforme pode sugerir.

Priorizar Popt, P95, Pfeas, Paceitavel, melhor valor, avaliações clássicas, shots e portas transpiladas em conjunto. Um critério de utilidade deve ser definido antes de ampliar a busca. Não elevar p apenas por expectativa de melhoria: o mixer factível ganhou pouco de p=2 para p=3, com muito mais portas.

Próximos testes que justificam investimento: preparação factível usando aritmética reversível/comparadores em vez de StatePreparation genérica; otimização progressiva por profundidade e maior orçamento de busca; reotimização com ruído e modelo calibrado; avaliação em instâncias novas, com custos, ameaças e criticidade diferentes. Não garantir que tais mudanças melhorarão esta ou outras instâncias. Corrigir amostras via heurística pode ser útil, mas os resultados precisam ser rotulados como pós-processamento clássico e comparados com a mesma heurística aplicada a amostras aleatórias.

As 84 execuções são explorações em uma mesma base ilustrativa, com poucas sementes; não são validação fora da amostra. O custo CVaR foi calculado com distribuições exatas, e nenhum experimento reotimizou parâmetros sob ruído. Sem evidência de vantagem quântica.

## Referências

- [Barkoutsos et al.: CVaR em otimização variacional](https://arxiv.org/abs/1907.04769).
- [Egger et al.: warm-starting quantum optimization](https://arxiv.org/abs/2009.10095).
- [Hadfield et al.: mixers para restrições](https://arxiv.org/abs/1709.03489).
- [Bärtschi e Eidenbenz: Grover mixers](https://arxiv.org/abs/2006.00354).

## Reproduzir

```bash
python -m unittest discover -v
python melhorias.py --ruido
python mixer_factivel.py --ruido
python analisar_resultados.py
```
