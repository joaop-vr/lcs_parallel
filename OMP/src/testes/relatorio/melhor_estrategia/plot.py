#!/usr/bin/env python3

import sys
import re
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def parse_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    pattern = re.compile(
        r'Configuração: threads=(\d+), sched=(\w+), chunk=(\d+)\.txt\s+'
        r'Média dos tempos totais:\s+([\d.]+) segundos',
        re.MULTILINE
    )

    dados = []
    for match in pattern.finditer(content):
        threads, sched, chunk, media = match.groups()
        dados.append({
            'threads': int(threads),
            'sched': sched,
            'chunk': int(chunk),
            'media': float(media)
        })

    return dados

def plotar_barras(dados):
    # Organiza os dados por número de threads
    categorias_sched_chunk = sorted(set(f"{d['sched']}_chunk{d['chunk']}" for d in dados))
    categorias_threads = sorted(set(d['threads'] for d in dados))

    # Mapeia tempos: {threads -> [tempos por config]}
    tempos_por_thread = defaultdict(lambda: [None]*len(categorias_sched_chunk))

    for d in dados:
        idx = categorias_sched_chunk.index(f"{d['sched']}_chunk{d['chunk']}")
        tempos_por_thread[d['threads']][idx] = d['media']

    x = np.arange(len(categorias_sched_chunk))  # posições no eixo X
    largura_barra = 0.18

    fig, ax = plt.subplots(figsize=(12, 6))

    cores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']  # até 4 cores para threads
    for i, thread in enumerate(categorias_threads):
        deslocamento = (i - len(categorias_threads)/2) * largura_barra + largura_barra/2
        tempos = tempos_por_thread[thread]
        ax.bar(x + deslocamento, tempos, width=largura_barra, label=f'{thread} threads', color=cores[i % len(cores)])

    ax.set_ylabel('Tempo médio (s)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(categorias_sched_chunk, rotation=45)
    ax.legend(title='Threads')
    ax.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Salvar o gráfico
    plt.savefig("grafico_barras_com_variancia.png", dpi=300)
    plt.show()

def main():
    if len(sys.argv) != 2:
        print(f'Uso: {sys.argv[0]} <arquivo_resultado.txt>')
        sys.exit(1)

    dados = parse_file(sys.argv[1])
    plotar_barras(dados)

if __name__ == '__main__':
    main()
