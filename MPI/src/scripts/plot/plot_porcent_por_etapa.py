#!/usr/bin/env python3
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Regex de extração de tempos
COLLECT_PATTERN = re.compile(r'Tempo de coleta da matriz,\s*([0-9.]+),')
TOTAL_PATTERN   = re.compile(r'Tempo total,\s*([0-9.]+),')

def extrair_coleta(path):
    text = open(path, 'r', encoding='utf-8').read()
    m_coll = COLLECT_PATTERN.search(text)
    m_tot  = TOTAL_PATTERN.search(text)
    t_coll = float(m_coll.group(1)) if m_coll else 0.0
    t_tot  = float(m_tot.group(1))  if m_tot  else t_coll
    pct    = 100.0 * t_coll / t_tot if t_tot > 0 else 0.0
    return pct

def coletar_dados(root_dir):
    configs, pct_collect = [], []
    pattern = re.compile(r'^[456]tasks$')        # só 4tasks, 5tasks ou 6tasks
    for sub in sorted(os.listdir(root_dir)):
        # filtra apenas os subdiretórios cujo nome bate com o regex
        if not pattern.match(sub):
            continue

        d = os.path.join(root_dir, sub)
        fp = os.path.join(d, 'metricas_output_teste_40000A_teste_40000B.out')
        if os.path.isdir(d) and os.path.isfile(fp):
            configs.append(sub)
            pct_collect.append(extrair_coleta(fp))
    return configs, np.array(pct_collect)

def plot_collect_only(configs, pct_collect, output_path):
    ind = np.arange(len(configs))
    width = 0.6

    plt.figure(figsize=(10,6))
    plt.bar(ind, pct_collect, width, color='C2')
    plt.ylabel('Porcentagem do tempo de coleta (%)')
    plt.title('Participação da Coleta no Tempo Total\nTeste 40000×40000')
    plt.xticks(ind, configs, rotation=45, ha='right')

    # anotações com valor percentual
    for i, pct in enumerate(pct_collect):
        plt.text(ind[i], pct/2, f'{pct:.1f}%', ha='center', va='center', color='white', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f'Gráfico salvo em: {output_path}')

def main():
    parser = argparse.ArgumentParser(
        description='Porcentagem de tempo de coleta (40000×40000)'
    )
    parser.add_argument('root_dir', help='Diretório com subdirs de resultados')
    args = parser.parse_args()

    if not os.path.isdir(args.root_dir):
        parser.error(f'"{args.root_dir}" não é um diretório válido.')

    configs, pct_collect = coletar_dados(args.root_dir)
    if not configs:
        print('Nenhum resultado encontrado em:', args.root_dir)
        return

    out_png = os.path.join(args.root_dir, 'percentual_coleta.png')
    plot_collect_only(configs, pct_collect, out_png)

if __name__ == '__main__':
    main()
