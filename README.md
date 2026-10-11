# Alocação de recursos de cibersegurança com QAOA

Implementação da formulação do artigo fornecido e dos experimentos do exercício OQA. Os PDFs são referências do problema; os dados numéricos em `instancia.json` são sintéticos e não foram extraídos de um levantamento real.

## Executar

```bash
git clone https://github.com/Lixipluv/CybersecQAOAMoney.git
cd CybersecQAOAMoney
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m unittest -v test_modelo.py
python projeto.py --ruido
```

Piloto rápido (cinco tamanhos, p=1/2/3 na maior instância e ruído na menor):

```bash
python projeto.py --ruido --repeticoes 1 --restarts 1 --maxiter 60 --shots 512 --saida resultados_piloto
```

O padrão usa três sementes e dois pontos iniciais por execução. `maxiter` é o limite de avaliações do COBYLA; a curva mostra avaliações, não iterações. Convergência e mensagem do otimizador ficam no JSON. O ruído usa os mesmos parâmetros otimizados no ideal, permitindo comparar o efeito da execução, sem reotimizar sob ruído.

## Formulação

Uma decisão binária por controle: x_i=1 significa selecionar o controle. O valor R_i implementa as equações 1–5 do artigo: eficácia multiplicada pela adesão humana; cobertura ponderada de táticas; prioridade CIS combinada com criticidade propagada e exposição de ativos; e baixo impacto operacional.

Maximizar sum(R_i x_i), sujeito a sum(c_i x_i) <= B. Custos e orçamento são inteiros na mesma unidade. Para dados monetários, escolha explicitamente uma unidade de discretização e avalie o erro introduzido.

A folga s=sum(2^m y_m), com M=ceil(log2(B+1)), transforma o orçamento em igualdade:

Q(x,y) = -sum(R_i x_i) + rho (sum(c_i x_i) + s - B)^2.

O padrão rho=1+sum(R_i) garante que uma violação inteira unitária não compense o ganho máximo possível. A representação matricial é triangular superior, com energia b^T Q b + rho B². O arquivo NPZ armazena esse offset e a constante Ising separadamente. Folgas maiores que B têm energia penalizada; não são soluções válidas de igualdade.

Com b_i=(1-z_i)/2, H=C+sum(h_i Z_i)+sum(J_ij Z_i Z_j). O circuito implementa H inicial, RZ(2 gamma h_i/escala), RZZ(2 gamma J_ij/escala) e RX(2 beta). A escala é max(abs(energias),1), preserva o ótimo e é registrada; gamma corresponde ao Hamiltoniano dividido por essa escala. A simulação NumPy é validada contra o statevector Qiskit.

## Experimentos e saídas

- Cinco tamanhos: 3–7 controles. O total de qubits inclui a folga, portanto difere do número de controles. Os subconjuntos são prefixos fixos da lista, sem aleatorização das instâncias.
- O orçamento de cada tamanho é floor(0.45 * custo de todos os controles), controlado por `fracao_orcamento`. O campo `orcamento` serve para uso direto das funções e testes; a varredura recalcula esse valor.
- Referência exata por enumeração de todos os estados. O simulador limita a 20 qubits; memória e tempo crescem exponencialmente.
- Profundidades p=1,2,3 na maior instância.
- Ruído opcional na menor instância: erros despolarizantes de 0.001 nas portas SX/X e 0.01 nas CX; erro de leitura de 0.02. São parâmetros sintéticos, sem calibração de dispositivo, sem decoerência temporal.
- Transpilação para conectividade linear, base RZ/SX/X/CX. Circuitos lógicos e transpilados são exportados em desenhos de texto; recursos após transpilação incluem medições na profundidade, mas excluem medições/barreiras da contagem de gates.
- Cada execução exporta modelo QUBO/Ising, parâmetros, energias, tempos, contagens, histogramas e curva de energia. O agregado contém JSON, CSV e seis gráficos de escalabilidade, com média e desvio padrão das repetições.

Popt agrega todas as codificações de portfólios ótimos do problema original, independentemente da folga. Pfeas mede orçamento válido no portfólio decodificado. E_melhor_medido mede a menor energia QUBO observada, que depende também da folga; a melhor solução factível por valor é registrada separadamente. Probabilidades exatas ideais e frequências de shots são identificadas separadamente. Não confundir energia esperada com energia do melhor bitstring.

Os mínimos de valor e cobertura são avaliados depois da otimização, como no artigo. O código registra `aceitabilidade_possivel`, `Paceitavel` e a aceitabilidade de cada solução relatada. Não impõe esses mínimos no QUBO. Valores de segurança aditivos podem contar benefícios redundantes; cobertura é a união das táticas, não a soma.

## Justificativa do Alinhamento com os Objetivos de Desenvolvimento Sustentável

ODS 16 — Paz, Justiça e Instituições Eficazes. O cenário ilustrativo protege portal público, dados de cidadãos e estações de atendimento. Selecionar controles dentro de um orçamento representa uma decisão institucional concreta. Maximizar o valor multicritério e avaliar cobertura busca apoiar proteção e continuidade dos serviços, confiança e tratamento responsável de dados. O indicador calculado é valor de segurança e cobertura ponderada, não redução de risco empiricamente medida nem impacto social comprovado.

Antes da entrega acadêmica, definir um critério de solução útil (por exemplo, valor próximo ao ótimo, probabilidade factível e aceitabilidade), executar uma exploração maior de tamanhos e orçamento, e discutir o limite observado. Sete controles é o limite desta base ilustrativa, não um limite demonstrado do QAOA. O piloto é uma verificação funcional; não demonstra vantagem quântica. A análise Q1–Q22 e conclusões sobre ruído, p e escalabilidade devem usar as saídas efetivamente obtidas.

Referências de implementação: https://qiskit.github.io/qiskit-aer/tutorials/3_building_noise_models.html e https://quantum.cloud.ibm.com/docs/en/guides/represent-quantum-computers . Os documentos originais permanecem nas localizações fornecidas pelo usuário.

## Melhorias experimentais e análise

A comparação de 84 execuções está em [ANALISE_MELHORIAS.md](ANALISE_MELHORIAS.md). O baseline original permanece como padrão. Os experimentos complementares implementam penalidade menor, CVaR, escala dos coeficientes, warm start heurístico e mixer factível de Grover. O ganho ideal do mixer é acompanhado de muito mais portas e piora sob o ruído testado; ele usa uma preparação genérica exponencial.

```bash
python -m unittest discover -v
python melhorias.py --ruido
python mixer_factivel.py --ruido
python analisar_resultados.py
python classico.py
```

`classico.py` resolve a mochila exatamente por programação dinâmica. `resultados_melhorias`, `resultados_mixer` e `resultados_analise` preservam parâmetros, distribuições, métricas e circuitos. As probabilidades ideais são exatas; as de ruído são frequências de shots. As análises mantêm o ótimo clássico fora da otimização variacional.

A função `modelo` também admite `rho > max(R)` para custos inteiros positivos e valores não negativos; a prova específica desta mochila está no relatório. O padrão `rho=1+sum(R)` foi preservado para reproduzir os resultados originais. O código de melhorias usa `rho=max(R)+0.05` nas variantes.

## Redução de decisões — caminho recomendado

A entrada `otimizar.py` usa **presolve certificado → QAOA no núcleo → reconstrução → verificação exata**. O orçamento utilizado é o do arquivo de dados, sem recalculá-lo automaticamente.

```bash
python otimizar.py --ruido
python otimizar.py --dados instancia.json --metodo rqaoa
python otimizar.py --dados instancia.json --metodo qiro
python otimizar.py --dados instancia.json --metodo agrupamento_qaoa --grupos 3
python experimentos_reducao.py --ruido
python auditar_reducao.py
```

Métodos disponíveis: `qaoa_completo`, `presolve_qaoa` (padrão), `agrupamento_qaoa`, `rqaoa`, `qiro`. `solucao_reduzida/resultado.json` contém a solução de exemplo; [relatório de redução](resultados_reducao/RELATORIO_REDUCAO.md) e [comparação](resultados_reducao/comparacao.png) documentam 150 execuções em dez cenários, com três sementes.

O presolve verifica desigualdades entre limites fracionários e portfólios factíveis usando aritmética racional. A redução preserva todos os ótimos da mochila valor/orçamento, inclusive empates. O núcleo mantém mapa de controles, offset do valor e orçamento residual; custos são simplificados pelo MDC. Se não restar nenhuma decisão, a solução é clássica, com zero qubits. O fato de o presolve preservar o ótimo não garante que a amostragem QAOA o encontre; o campo `otima_verificada` é uma comparação posterior com DP racional no modelo original.

RQAOA e QIRO usam correlações estimadas por shots ideais para reduções heurísticas, com cauda clássica exata. QIRO é uma adaptação para mochila com guarda de orçamento e backtracking limitado, sem garantia global. O agrupamento escolhe representantes por valor/custo e pode eliminar controles necessários ao ótimo. Esses três métodos são experimentais: as perdas observadas ficam explícitas no relatório. Não se aplica dominância simples entre controles que podem ser selecionados juntos.

Validação: 22 testes, incluindo 120 mochilas aleatórias com empates, verificação de limites, energia/reconstrução Ising, orçamento zero, certificados adulterados e ruído no núcleo. A auditoria das 150 execuções confirmou 138 fixações certificadas e 30 reconstruções RQAOA. Nessa amostra, o pipeline certificado obteve o ótimo em 30/30 execuções, sendo 21 resolvidas classicamente pelo presolve. Dados sintéticos; sem alegação de vantagem quântica.

## Gráficos da evolução

A [galeria comparativa](graficos_evolucao/README.md) reúne cinco gráficos em PNG e SVG: evolução histórica, antes/depois sob ruído, recursos por tamanho, sucesso das reduções e contribuição clássica. Os dados usados estão disponíveis em CSV e JSON.

```bash
python graficos_evolucao.py
```

Os gráficos distinguem probabilidades por medição, melhor solução encontrada e sucesso por execução. A comparação histórica identifica diferenças de configuração; o antes/depois do presolve usa configurações controladas.
