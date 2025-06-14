#!/bin/bash

# Diretório para armazenar os resultados (cria se não existir)
OUTPUT_DIR="resultados"
mkdir -p "$OUTPUT_DIR"

# Tamanhos possíveis (mais granular, mas ainda respeitando os limites)
INPUT_SIZES=(1000 2500 5000 7500 10000 25000 50000 75000 100000)

# Threads, escalonamentos e chunks
THREADS=(1 2 4 8)
SCHEDULES=("static" "dynamic" "guided")
CHUNKS=(64 128)

# Loop principal (combinando todos os pares diferentes de tamanhos A e B)
for sizeA in "${INPUT_SIZES[@]}"; do
    for sizeB in "${INPUT_SIZES[@]}"; do

        # pula se os tamanhos forem iguais
        if [[ "$sizeA" -eq "$sizeB" ]]; then
            continue
        fi

        # Gera arquivos de entrada distintos
        ./generate_input "$sizeA" "$sizeB" "teste_${sizeA}_${sizeB}"

        for threads in "${THREADS[@]}"; do
            for chunk in "${CHUNKS[@]}"; do
                for schedule in "${SCHEDULES[@]}"; do
                    echo "[TESTE] A=$sizeA | B=$sizeB | threads=$threads | schedule=$schedule | chunk=$chunk"

                    OUTPUT_FILE="${OUTPUT_DIR}/A_${sizeA}_B_${sizeB}_threads_${threads}_sched_${schedule}_chunk_${chunk}.txt"

                    for run in {1..21}; do
                        echo "--- Execução $run ---" >> "$OUTPUT_FILE"
                        ./lcs "teste_${sizeA}_${sizeB}A.in" "teste_${sizeA}_${sizeB}B.in" "$threads" "$chunk" "$schedule" >> "$OUTPUT_FILE"
                        echo "" >> "$OUTPUT_FILE"
                    done
                done
            done
        done

        # rm teste_${sizeA}_${sizeB}*.in  # descomente se quiser limpar após
    done
done

echo "Todos os testes concluídos! Resultados salvos em: $OUTPUT_DIR/"
