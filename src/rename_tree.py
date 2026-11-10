#!/usr/bin/env python3
"""Переименовывает листья Newick-дерева по таблице genomes.tsv.

Использование:
    python rename_tree.py <папка>

В <папка>/trees берётся файл дерева без "full" в названии; метки заменяются
на значения из последней колонки <папка>/genomes.tsv. Результат сохраняется
в <папка>/"trees renamed"/<имя исходного файла>.
"""
import re
import sys
from pathlib import Path
import argparse
import pandas as pd

# лист: метка после "(" или ",", за которой идёт ":" (длина ветви) или закрывающая скобка/запятая
LEAF_RE = re.compile(r"(?P<pre>[(,]\s*)(?P<label>'[^']*(?:''[^']*)*'|[^():,;\s]+)(?P<post>\s*(?=[:,)]))")

def quote(name: str) -> str:
    """Убирает запятые, заменяет пробелы на "_", остальные спецсимволы Newick берёт в кавычки."""
    name = re.sub(r"\s+", "_", name.replace(",", "").strip())
    if re.search(r"[\s():;,\[\]']", name):
        return "'" + name.replace("'", "''") + "'"
    return name

def parse_args():
    parser = argparse.ArgumentParser(prog = "rename_tree", description="Изменяет лейблы на деревьях")
    parser.add_argument("indir", help = "Директория пангенома", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    folder = args.indir
    tree_path = folder/"trees"/"nj-global-tree.tre"
    text = tree_path.read_text(encoding="utf-8")

    genomes = pd.read_csv(folder / "genomes.tsv", sep="\t", header=None, names = ["ac", "abbrev", "seq_type", "c/l", "definition"])
    mapping = dict(zip(genomes["abbrev"], genomes["definition"]))

    def repl(m: re.Match) -> str:
        label = m.group("label").strip("'")
        return m.group("pre") + quote(mapping[label]) + m.group("post")

    renamed = LEAF_RE.sub(repl, text)

    out_dir = folder / "trees renamed"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / tree_path.name
    out_path.write_text(renamed, encoding="utf-8")
    print(f"Сохранено: {out_path}")


if __name__ == "__main__":
    main()
