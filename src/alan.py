import argparse
import os
import pandas as pd
import copy as cp
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(prog="AlAn", description="AlAn, version 2.0.0")
    parser.add_argument("indir", nargs="?", default=os.getcwd(),
                         help="input directory (default - current)")
    parser.add_argument("outdir", nargs="?", default=os.getcwd(),
                         help="output directory (default - current)")
    parser.add_argument("-c", dest="out_excel", action="store_false", default=True,
                         help="output files type (default - excel, otherwise - csv)")
    return parser.parse_args()

def write_table(table: pd.DataFrame, name: str, out: str, out_excel: bool, *, index: bool = True) -> None:
    if out_excel:
        table.to_excel(f"{out}/{name}.xlsx", index = index)
    else:
        table.to_csv(f"{out}/{name}.csv", index = index)

def group_table_f(group_table_in: pd.DataFrame, group_id: int, outdir: str, mut_path: str, out_excel: bool):
    result = {"gene": [], "ori": [], "block_start": [], "block_stop": [], "gene_start": [], "gene_stop": [], "len": [], "start_rf": [], "stop_rf": [], "is_start": [], "is_stop": [], "past_group": [], "next_group": []}
    for index_counter in range(len(group_table_in)):
        result["gene"].append(group_table_in.iloc[index_counter, 3])
        result["ori"].append(group_table_in.iloc[index_counter, 7])
        result["block_start"].append(group_table_in.iloc[index_counter, 5] + 1)
        result["block_stop"].append(group_table_in.iloc[index_counter, 6] + 1)
        result["gene_start"].append(group_table_in.iloc[index_counter, 8] + 1)
        result["gene_stop"].append(group_table_in.iloc[index_counter, 9] + 1)
        result["len"].append(group_table_in.iloc[index_counter, 6] - group_table_in.iloc[index_counter, 5] + 1)
        result["is_start"].append(group_table_in.iloc[index_counter, 8] == 0)
        result["is_stop"].append(group_table_in.iloc[index_counter, 10])
        try:
            with open(f"{mut_path}/{group_table_in.iloc[0, 4]}.tsv") as mutfile:
                succ = False
                while not succ:
                    line = mutfile.readline()
                    line_split = line.split("\t")
                    if line == "":
                        succ = True
                    elif line_split[1] != ".":
                        seq_split = line_split[1].split("_")
                        if group_table_in.iloc[index_counter, 1] <= group_table_in.iloc[index_counter, 2]:
                            if (seq_split[0] == group_table_in.iloc[index_counter, 0]) and (group_table_in.iloc[index_counter, 1] >= int(seq_split[1])) and \
                            (group_table_in.iloc[index_counter, 2] <= int(seq_split[2])):
                                succ = True
                        else:
                            if (seq_split[0] == group_table_in.iloc[index_counter, 0]) and (group_table_in.iloc[index_counter, 1] <= int(seq_split[1])) and \
                            (group_table_in.iloc[index_counter, 2] >= int(seq_split[2])):
                                succ = True
                split_sum = 0
                if line != "":
                    coord = int(line_split[2])
                    if line_split[3].strip() == "-":
                        split_sum += 1
                    elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                        split_sum += int(line_split[3]) - coord + 1
                    line = mutfile.readline()
                    line_split = line.split("\t")
                    if line != "":
                        coord = int(line_split[2])
                while (line != "") and (coord < group_table_in.iloc[index_counter, 5]) and (line_split[1] == "."):
                    if line_split[3].strip() == "-":
                        split_sum += 1
                    elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                        split_sum += int(line_split[3]) - coord + 1
                    line = mutfile.readline()
                    line_split = line.split("\t")
                    if line != "":
                        coord = int(line_split[2])
                if group_table_in.iloc[index_counter, 7] == 1:
                    result["start_rf"].append((group_table_in.iloc[index_counter, 5] + 1 - split_sum)%3 + 1)
                else:
                    result["start_rf"].append(-((group_table_in.iloc[index_counter, 5] + 1 - split_sum)%3 + 1))
                while (line != "") and (coord < group_table_in.iloc[index_counter, 6]) and (line_split[1] == "."):
                    if line_split[3].strip() == "-":
                        split_sum += 1
                    elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                        split_sum += int(line_split[3]) - coord + 1
                    line = mutfile.readline()
                    line_split = line.split("\t")
                    if line != "":
                        coord = int(line_split[2])
                if group_table_in.iloc[index_counter, 7] == 1:
                    result["stop_rf"].append((group_table_in.iloc[index_counter, 6] - 1 - split_sum)%3 + 1)
                else:
                    result["stop_rf"].append(-((group_table_in.iloc[index_counter, 6] - 1 - split_sum)%3 + 1))
        except FileNotFoundError:
            if group_table_in.iloc[index_counter, 7] == 1:
                result["start_rf"].append((group_table_in.iloc[index_counter, 5] + 1)%3 + 1)
                result["stop_rf"].append((group_table_in.iloc[index_counter, 6] - 1)%3 + 1)
            else:
                result["start_rf"].append(-((group_table_in.iloc[index_counter, 5] + 1)%3 + 1))
                result["stop_rf"].append(-((group_table_in.iloc[index_counter, 6] - 1)%3 + 1))
        result["past_group"].append(np.nan)
        result["next_group"].append(np.nan)
    group_table_out = pd.DataFrame(result)
    linked_table = group_table_out[~group_table_out["is_start"] | ~group_table_out["is_stop"]]
    linked_table.insert(13, "group", [group_id]*len(linked_table))
    write_table(group_table_out, group_id, outdir, out_excel)
    return(group_table_out["is_start"].all(), group_table_out["is_stop"].all(), linked_table)

def final_statistics(table: pd.DataFrame, outdir: str, out_excel: bool) -> None:
    counts = table["group_type"].value_counts().sort_index().reset_index()
    counts.columns = ["group_type", "count"]
    write_table(counts, "final_statistics", outdir, out_excel, index = False)

def name_change(groups: pd.DataFrame, indir: str) -> pd.DataFrame:
    genomes = pd.read_csv(f"{indir}/genomes.tsv", sep = "\t", header = None, names = ["ac", "abbrev", "seq_type", "c/l", "definition"])
    genomes["merged_cols"] = genomes["abbrev"] + "&" + genomes["seq_type"] +  "&" + genomes["c/l"]
    mapping = dict(zip(genomes["merged_cols"], genomes["definition"]))
    groups.columns = [mapping.get(col, col) for col in groups.columns]
    return groups

def main():
    args = parse_args()
    indir = args.indir
    outdir = args.outdir
    out_excel = args.out_excel

    os.makedirs(outdir, exist_ok=True)

    table = pd.read_csv(f"{indir}/genes/partition-ungrouped.tsv", sep = "\t")
    table = table.iloc[:, :-1]
    table["gene_end"] = [False]*len(table)
    for gene, loc_table in table.groupby("gene"):
        loc_table_sort = loc_table.sort_values("gene_block_stop", ascending = False)
        table.loc[loc_table_sort.iloc[0].name, "gene_end"] = True
    try:
        os.mkdir(f"{outdir}/mut_split")
    except FileExistsError:
        pass
    with open(f"{indir}/mutations/mut.tsv") as infile:
        line = infile.readline()
        line = infile.readline()
        line_split = line.split("\t")
        while line != "":
            with open(f"{outdir}/mut_split/{line_split[0]}.tsv", "w") as outfile:
                outfile.write(line)
                line = infile.readline()
                line_split = line.split("\t")
                while (line != "") and (line_split[0] == "."):
                    outfile.write(line)
                    line = infile.readline()
                    line_split = line.split("\t")

    result = {"block": [], "block_type": [], "seq_num": [], "block_len": [], "ori": [], "start": [], "stop": [], "len": [], "core_start": [], "core_stop": [], "core_len": [], "genes_num": [], "is_true_start": [], "is_true_stop": [], "core_%": [], "group_type": []}
    seqs = table["sequence"].unique()
    for seq in seqs:
        result[seq] = []
    linked_table = pd.DataFrame({"gene": [], "ori": [], "block_start": [], "block_stop": [], "gene_start": [], "gene_stop": [], "len": [], "start_rf": [], "stop_rf": [], "is_start": [], "is_stop": [], "past_group": [], "next_group": [], "group": []})
    for block_ori, loc_table in table.groupby(["npg_block", "npg_block_ori"]):
        starts = sorted(loc_table["npg_block_min"])
        stops = sorted(loc_table["npg_block_max"])
        starts_stops = [[0, True]]*(len(starts) + len(stops))
        starts_counter = 0
        stops_counter = 0
        while starts_counter < len(starts):
            if starts[starts_counter] <= stops[stops_counter]:
                starts_stops[starts_counter + stops_counter] = [starts[starts_counter], True]
                starts_counter += 1
            else:
                starts_stops[starts_counter + stops_counter] = [stops[stops_counter], False]
                stops_counter += 1
        for stops_final_counter in range(stops_counter, len(stops)):
            starts_stops[starts_counter + stops_final_counter] = [stops[stops_final_counter], False]
        last_coord = -1
        for starts_stops_counter in range(len(starts_stops)):
            if starts_stops[starts_stops_counter][1]:
                last_coord = starts_stops[starts_stops_counter][0]
            elif last_coord != -1:
                core_start = last_coord
                core_stop = starts_stops[starts_stops_counter][0]
                group_table = loc_table[(loc_table["npg_block_min"] <= core_start) & (loc_table["npg_block_max"] >= core_stop)]
                group_start = min(group_table["npg_block_min"])
                group_stop = max(group_table["npg_block_max"])
                block_split = block_ori[0].replace("n", "x").split("x")
                result["block"].append(block_ori[0])
                result["block_type"].append(block_split[0][0])
                result["seq_num"].append(int(block_split[0][1:]))
                result["block_len"].append(int(block_split[1]))
                result["ori"].append(int(block_ori[1]))
                result["start"].append(group_start + 1)
                result["stop"].append(group_stop + 1)
                result["len"].append(group_stop - group_start + 1)
                result["core_start"].append(core_start + 1)
                result["core_stop"].append(core_stop + 1)
                result["core_len"].append(core_stop - core_start + 1)
                result["genes_num"].append(len(group_table))
                seqs_counter = 0
                for seq in seqs:
                    if seq in group_table["sequence"].unique():
                        result[seq].append(True)
                        seqs_counter += 1
                    else:
                        result[seq].append(False)
                last_coord = -1
                try:
                    os.mkdir(f"{outdir}/groups")
                except FileExistsError:
                    pass
                is_start, is_stop, new_linked_table = group_table_f(group_table, len(result["block"]) - 1, f"{outdir}/groups", f"{outdir}/mut_split", out_excel)
                try:
                    linked_table = linked_table.merge(new_linked_table, "outer")
                except Exception:
                    pass
                result["is_true_start"].append(is_start)
                result["is_true_stop"].append(is_stop)
                core_percent = round((core_stop - core_start + 1)/(group_stop - group_start + 1), 4)*100
                result["core_%"].append(core_percent)
                mistake_type = ""
                if ((result["ori"][-1] == 1) and (result["start"][-1] != result["core_start"][-1])) or ((result["ori"][-1] == -1) and (result["stop"][-1] != result["core_stop"][-1])):
                    mistake_type += "/wrong start"
                if (result["genes_num"][-1] < result["seq_num"][-1]) and is_start and is_stop:
                    mistake_type += "/annotation loss"
                if result["genes_num"][-1] > result["seq_num"][-1]:
                    mistake_type += "/crossing genes"
                if mistake_type == "":
                    mistake_type = "no mistake"
                result["group_type"].append(mistake_type)
    groups = pd.DataFrame(result)
    final_statistics(groups, outdir, out_excel)
    write_table(linked_table, "linked_table", outdir, out_excel)
    linked_groups = groups[groups["group_type"] == 3]
    for line_number in range(len(linked_groups)):
        if out_excel:
            group = pd.read_excel(f"{outdir}/groups/{linked_groups.iloc[line_number, :].name}.xlsx", index_col = 0)
        else:
            group = pd.read_csv(f"{outdir}/groups/{linked_groups.iloc[line_number, :].name}.csv")
        for gene_number in range(len(group)):
            if not group.iloc[gene_number, 9]:
                try:
                    group.iloc[gene_number, 11] = linked_table[(linked_table["gene"] == group.iloc[gene_number, 0]) & ((linked_table["gene_stop"] + 1) == group.iloc[gene_number, 4])].iloc[0, 13]
                except Exception:
                    pass
            if not group.iloc[gene_number, 10]:
                try:
                    group.iloc[gene_number, 12] = linked_table[(linked_table["gene"] == group.iloc[gene_number, 0]) & ((linked_table["gene_start"] - 1) == group.iloc[gene_number, 5])].iloc[0, 13]
                except Exception:
                    pass
        write_table(group, linked_groups.iloc[line_number, :].name, f"{outdir}/groups/", out_excel)
    groups = name_change(groups, indir)
    write_table(groups, "groups", outdir, out_excel)

if __name__ == "__main__":
    main()