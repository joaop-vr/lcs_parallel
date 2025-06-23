import os
import re
import argparse
import numpy as np

def process_file(path):
    """
    Lê o arquivo em 'path', extrai métricas numéricas e retorna um dicionário
    com listas de valores para cada métrica.
    Suporta formatos com ':' ou '=' e valores com unidades (ex: 'segundos').
    """
    metrics = {}
    # Regex para capturar 'nome_metrica: valor' ou 'nome_metrica = valor'
    pattern = re.compile(r"^\s*([\wáéíóúãõâêîôûçÁÉÍÓÚÃÕÂÊÎÔÛÇ\(\)\- ]+?)\s*[:=]\s*([0-9]+\.?[0-9]*)")
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            # ignora cabeçalhos de timestamp
            if line.strip().startswith('['):
                continue
            m = pattern.match(line)
            if m:
                name, val = m.groups()
                try:
                    num = float(val)
                except ValueError:
                    continue
                key = name.strip()
                # remover sufixos irrelevantes
                key = key.replace('segundos', '').strip()
                metrics.setdefault(key, []).append(num)
    return metrics


def compute_stats(metrics):
    """
    Recebe dicionário de listas e retorna um dicionário
    com (média, desvio_padrão_relativo) para cada métrica.
    """
    stats = {}
    for name, vals in metrics.items():
        arr = np.array(vals)
        mean = arr.mean()
        std = arr.std(ddof=0)
        rel_std = std / mean if mean != 0 else float('nan')
        stats[name] = (mean, rel_std)
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Calcular média e desvio padrão relativo de métricas em vários arquivos.')
    parser.add_argument(
        'directory',
        help='Diretório contendo arquivos a serem processados')
    args = parser.parse_args()

    for root, _, files in os.walk(args.directory):
        for fname in files:
            path = os.path.join(root, fname)
            metrics = process_file(path)
            if not metrics:
                continue  # pula arquivos sem métricas válidas

            stats = compute_stats(metrics)
            out_name = f"metricas_{fname}"
            out_path = os.path.join(root, out_name)
            with open(out_path, 'w', encoding='utf-8') as out:
                out.write(f"Estatísticas para {fname}\n")
                out.write("Métrica, Média, Desvio Padrão Relativo\n")
                for name, (mean, rel) in stats.items():
                    out.write(f"{name}, {mean:.6f}, {rel:.6f}\n")
            print(f"Gerado: {out_path}")

if __name__ == '__main__':
    main()
