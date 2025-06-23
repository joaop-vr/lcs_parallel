#!/usr/bin/env python3
import os
import re
import argparse
import matplotlib.pyplot as plt

TEMPO_TOTAL_PATTERN = re.compile(r'Tempo total,\s*([0-9.]+),\s*([0-9.]+)')

def coletar_dados(root_dir):
    dados = {}
    for subdir in sorted(os.listdir(root_dir)):
        path_sub = os.path.join(root_dir, subdir)
        if not os.path.isdir(path_sub):
            continue

        xs, ys, errs = [], [], []
        for fname in sorted(os.listdir(path_sub)):
            if not fname.endswith('.out'):
                continue

            m = re.search(r'teste_(\d+)A', fname)
            if not m:
                continue
            tamanho = int(m.group(1))

            with open(os.path.join(path_sub, fname), 'r', encoding='utf-8') as f:
                content = f.read()
            m2 = TEMPO_TOTAL_PATTERN.search(content)
            if not m2:
                continue
            tempo = float(m2.group(1))
            desvio = float(m2.group(2))

            xs.append(tamanho)
            ys.append(tempo)
            errs.append(desvio)

        if xs:
            sorted_data = sorted(zip(xs, ys, errs), key=lambda t: t[0])
            xs, ys, errs = zip(*sorted_data)
            dados[subdir] = {'x': xs, 'y': ys, 'yerr': errs}
    return dados

def plotar(dados, output_path):
    plt.figure(figsize=(8, 6))
    for subdir, dt in dados.items():
        plt.errorbar(
            dt['x'], dt['y'], yerr=dt['yerr'],
            marker='o', linestyle='-',
            label=subdir
        )
    plt.xlabel('Tamanho do teste (n)')
    plt.ylabel('Tempo total (s)')
    plt.title('Tempo total vs. Tamanho do teste')
    plt.legend(title='Número Processos')
    plt.grid(True)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300)
    print(f'✅ Gráfico salvo em: {output_path}')

def main():
    parser = argparse.ArgumentParser(
        description='Lê arquivos de métricas e salva gráfico do tempo total com desvio padrão.'
    )
    parser.add_argument(
        'root_dir',
        help='Caminho para o diretório contendo os subdiretórios (1task, 2tasks, …)'
    )
    args = parser.parse_args()

    if not os.path.isdir(args.root_dir):
        parser.error(f'O caminho "{args.root_dir}" não é um diretório válido.')

    dados = coletar_dados(args.root_dir)
    if not dados:
        print('Nenhum dado válido encontrado em:', args.root_dir)
        return

    output_path = os.path.join(args.root_dir, 'grafico_metricas.png')
    plotar(dados, output_path)

if __name__ == '__main__':
    main()
