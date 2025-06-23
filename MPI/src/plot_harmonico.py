#!/usr/bin/env python3
import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Regex de extração de tempos
P_PATTERN       = re.compile(r'Tempo de construção de P,\s*([0-9.]+),')
W_PATTERN       = re.compile(r'Tempo de cálculo DP \(wavefront\),\s*([0-9.]+),')
INIT_PATTERN    = re.compile(r'Tempo de inicialização \(dados\),\s*([0-9.]+),')
COLLECT_PATTERN = re.compile(r'Tempo de coleta da matriz,\s*([0-9.]+),')
TOTAL_PATTERN   = re.compile(r'Tempo total,\s*([0-9.]+),')

def extrair_tempos(path):
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()
    def busca(pat):
        m = pat.search(text)
        return float(m.group(1)) if m else 0.0
    t_p    = busca(P_PATTERN)
    t_w    = busca(W_PATTERN)
    t_init = busca(INIT_PATTERN)
    t_coll = busca(COLLECT_PATTERN)
    t_tot  = busca(TOTAL_PATTERN) or (t_p + t_w + t_init + t_coll)
    return t_p, t_w, t_init + t_coll, t_tot

def coletar_dados(root_dir):
    configs, comp, comm = [], [], []
    for subdir in sorted(os.listdir(root_dir)):
        path_sub = os.path.join(root_dir, subdir)
        if not os.path.isdir(path_sub):
            continue
        fp = os.path.join(path_sub, 'metricas_output_teste_40000A_teste_40000B.out')
        if not os.path.isfile(fp):
            continue
        t_p, t_w, t_c, t_tot = extrair_tempos(fp)
        configs.append(subdir)
        comp.append(t_p + t_w)
        comm.append(t_c)
    return configs, np.array(comp), np.array(comm)

def calc_score_harmonico(comp, comm):
    # evita divisão por zero
    return 2 * comp * comm / (comp + comm + 1e-16)

def plot_scatter(configs, comp, comm, scores, best_idx, second_idx, output_path):
    plt.figure(figsize=(8,6))
    for i, cfg in enumerate(configs):
        if i == best_idx or i == second_idx:
            size = 150 if i == best_idx else 100
            edge = 'black'
            alpha = 1.0
            label = f"{cfg} (score={scores[i]:.4f})"
            label = f"{cfg} {'(Best)' if i==best_idx else '(2nd best)'} score={scores[i]:.4f}"
            plt.scatter(comp[i], comm[i], s=size, edgecolor=edge, alpha=alpha, label=label)
        else:
            plt.scatter(comp[i], comm[i], s=60, alpha=0.6, label=cfg)
    plt.xlabel('Tempo de Computação (s)')
    plt.ylabel('Tempo de Comunicação (s)')
    plt.title('Trade-off Computação vs Comunicação (n=40000×40000)')
    plt.legend(loc='best', fontsize='small', ncol=2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f'✅ Gráfico salvo em: {output_path}')
    print(f'🏆 Melhor configuração: {configs[best_idx]} (score harmônico = {scores[best_idx]:.6f})')
    print(f'🥈 2ª melhor configuração: {configs[second_idx]} (score harmônico = {scores[second_idx]:.6f})')

def main():
    parser = argparse.ArgumentParser(
        description='Avalia trade-off e score harmônico para teste 40000×40000'
    )
    parser.add_argument('root_dir', help='Diretório com subdirs de tasks')
    args = parser.parse_args()
    
    if not os.path.isdir(args.root_dir):
        parser.error(f'"{args.root_dir}" não é um diretório válido.')

    configs, comp, comm = coletar_dados(args.root_dir)
    if not configs:
        print('Nenhum teste 40000×40000 encontrado em:', args.root_dir)
        return

    scores = calc_score_harmonico(comp, comm)
    sorted_idx = np.argsort(scores)
    best_idx = sorted_idx[0]
    second_idx = sorted_idx[1] if len(sorted_idx) > 1 else None

    output_png = os.path.join(args.root_dir, 'scatter_score_harmonico.png')
    plot_scatter(configs, comp, comm, scores, best_idx, second_idx, output_png)

if __name__ == '__main__':
    main()
