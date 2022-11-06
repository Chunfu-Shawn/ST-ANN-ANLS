import sys
import argparse
import numpy as np
import utils
import pandas as pd


def rename_gene(input_names,
                input_species="mouse"):
    """transfer mouse gene names to human ortholog gene names

    Parameters
    ---------
    input_names: string or np.ndarray
        input gene names
    input_species: string
        transfer mouse to human ortholog genes or rather

    Returns
    ---------
    genes: np.ndarray
        tranformed gene name or gene ensembl id

    """
    utils.warn(" --- Load Data... ---")
    if type(input_names) == "str":
        genes = input_names.split(",")
    else:
        genes = input_names
    print(genes)
    # input_data = pd.read_csv(input_file, sep=',')
    ref_mouse2human = "/home/user/data2/rbase/spatial_annotate_scripts/src/orthologs_mouse2human.csv"
    ref_human2mouse = "/home/user/data2/rbase/spatial_annotate_scripts/src/orthologs_human2mouse.csv"
    ref_id = {}
    ref_syb = {}
    utils.warn(" --- Create Reference Orthologues between Human and Mouse ---")
    if input_species == "mouse":
        for line in open(ref_mouse2human, "r"):
            # gene: mouse_id mouse_name human_id human_name
            line = line.strip()
            line = line.split(",")
            m_gene_id = line[0]
            m_gene_syb = line[1]
            h_gene_id = line[2]
            h_gene_syb = line[3]
            ref_id[m_gene_id] = h_gene_id
            ref_syb[m_gene_syb] = h_gene_syb
    elif input_species == "human":
        for line in open(ref_human2mouse, "r"):
            # gene: human_id human_name mouse_id mouse_name
            line = line.split(",")
            h_gene_id = line[0]
            h_gene_syb = line[1]
            m_gene_id = line[2]
            m_gene_syb = line[3]
            ref_id[h_gene_id] = m_gene_id
            ref_syb[h_gene_syb] = m_gene_syb
    utils.warn(" --- Begin Transferring ---")
    for i in range(0, len(genes)):
        if genes[i] in ref_id.keys() and ref_id[genes[i]] != '':
            genes[i] = ref_id[genes[i]]
            print(genes[i])
        elif genes[i] in ref_syb.keys() and ref_syb[genes[i]] != '':
            genes[i] = ref_syb[genes[i]]
            print(genes[i])
    return genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Rename mouse gene name to human orthologs gene name or reversely')
    input_str = sys.stdin.read()
    parser.add_argument('--input_species', '-s', required=True, help='input gene name is mouse or human')
    args = parser.parse_args()
    print(rename_gene(input_names=input_str, **vars(args)))
