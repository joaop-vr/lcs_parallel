## Plano de Testes para Avaliação do LCS Paralelo

Este documento descreve detalhadamente o plano de testes para atender aos requisitos do trabalho de paralelização do algoritmo LCS.

---

### 1. Ambiente Experimental

* **Máquina:** Dell PowerEdge R740
* **Processador:** Intel Xeon E5-2670 v3 (12 cores físicos, 24 HT)
* **Sistema Operacional:** Ubuntu 22.04 LTS, kernel 5.15.x
* **Compilador:** GCC 12.2.0
* **Flags de Compilação:** `-O3 -march=native -fopenmp`
* **Configurações Extras:**

  * Desativar turbo-boost / DVFS:

    ```bash
    echo performance > /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
    ```
  * Pinagem de threads (affinity):

    ```bash
    export OMP_PROC_BIND=spread
    export OMP_PLACES=cores
    ```

---

### 2. Metodologia Geral

1. **Corretude**: Verificar se a saída paralela bate com a sequencial para diversas entradas (aleatórias e casos específicos).
2. **Repetições**: Executar cada combinação de parâmetros (threads, tamanhos) 20 vezes para obter média e desvio-padrão.
3. **Medição de tempo**: Utilizar `omp_get_wtime()` nos seguintes trechos:

   * **T sequencial puro**: trecho inicial antes da primeira diretiva.
   * **T paralelizável**: tempo total menos o tempo sequencial.
4. **Coleta de dados**: Gerar logs no formato CSV:

   ```csv
   modo,tamanho,threads,run,tempo_total,tempo_seq,tempo_par
   seq,10000,1,1,12.345,12.345,0.000
   par,10000,4,1,3.210,0.500,2.710
   ```

---

### 3. Teste de Escalabilidade Forte

* **Tamanho fixo**: escolher N tal que o tempo sequencial seja ≥ 10 segundos (ex.: N = 200000).
* **Variar número de threads**: 1, 2, 4, 8, 12, 24.
* **Cálculos**:

  * **Speedup**: $S_p = T(1) / T(p)$
  * **Eficiência**: $E_p = S_p / p$

---

### 4. Teste de Escalabilidade Fraca

* **Aumentar N proporcionalmente** ao número de threads: se N₁ = 100000 para 1 thread, usar $N_p = 100000 \times p$.
* **Variar p**: 1, 2, 4, 8, 12, 24.
* **Speedup fraco**: $S_p^{(fraco)} = T(1, N₁) / T(p, N_p)$.

---

### 5. Lei de Amdahl

1. Estimar fração sequencial $f = \frac{T_{seq}}{T_{total}}$ em 1 thread.
2. Calcular speedup teórico:
   $S_p = \frac{1}{f + \frac{1-f}{p}}$
3. Tabela de speedup teórico para p = 2, 4, 8 e ∞:

| p | Sₚ teórico      |
| - | --------------- |
| 2 | $1/(f+(1-f)/2)$ |
| 4 | $1/(f+(1-f)/4)$ |
| 8 | $1/(f+(1-f)/8)$ |
| ∞ | $1/f$           |

---

### 6. Tabelas no Relatório

* **Escalabilidade Forte**: tabela de speedup e eficiência versus p para N fixo.
* **Escalabilidade Fraca**: tabela de speedup fraco versus p.

---

### 7. Análise dos Resultados

* Comparar speedup obtido com speedup linear e com Lei de Amdahl.
* Discutir overheads (sincronização, carga desigual nas diagonais).
* Verificar ruídos nos dados (desvio-padrão alto) e possíveis causas.

---

### 8. Escalabilidade

* **Forte**: avaliar se o speedup cresce de forma sublinear próxima a p.
* **Fraca**: avaliar se a eficiência se mantém estável conforme p aumenta.

---

### 9. Observações Finais

* Garantir que ambos os códigos (seq. e paralelo) gerem saídas idênticas.
* Executar em ambiente dedicado para reduzir variabilidade.
* Documentar versão de software e hardware usados.

---

*Este plano cobre todos os itens do enunciado e fornece base estatística sólida para o relatório.*
