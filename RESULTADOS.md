# Execução validada

Foram executados 24 experimentos com três sementes (42–44), 2048 shots, dois pontos iniciais e até 150 avaliações COBYLA por ponto. Dados sintéticos. Cinco testes passaram, incluindo equivalência integral QUBO/Ising, circuito NumPy/Qiskit, leitura após transpilação e exportação dos resultados.

| Controles | p | Ambiente | Qubits | Popt média | Pfeas média |
|---|---|---|---|---|---|
| 3 | 1 | ideal | 5 | 0.2052 | 0.7928 |
| 3 | 1 | ruido_sintetico | 5 | 0.1929 | 0.7357 |
| 4 | 1 | ideal | 7 | 0.0672 | 0.7393 |
| 5 | 1 | ideal | 8 | 0.0400 | 0.7803 |
| 6 | 1 | ideal | 9 | 0.0146 | 0.7896 |
| 7 | 1 | ideal | 11 | 0.0101 | 0.7896 |
| 7 | 2 | ideal | 11 | 0.0065 | 0.8136 |
| 7 | 3 | ideal | 11 | 0.0112 | 0.8221 |

A frequência ótima diminuiu com o tamanho nesta família. Em sete controles, aumentar p não trouxe melhoria consistente de Popt; a factibilidade aumentou em média. O ruído sintético reduziu a factibilidade na instância de três controles. Essas observações são específicas destes dados, sementes e orçamento de otimização.

Os dados não determinam o maior tamanho útil: é necessário definir o critério de utilidade e ampliar a família. A conexão com ODS 16 é operacional e ilustrativa; não há medição de impacto social ou vantagem quântica.

Para reproduzir com as mesmas versões, use `pip install -r requirements-lock.txt`.
