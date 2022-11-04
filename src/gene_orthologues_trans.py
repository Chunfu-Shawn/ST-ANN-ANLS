import sys
import os
import csv
import argparse
import utils
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(sys.argv[0])))


def rename_gene(input_names, ref="/home/user/data2/rbase/spatial_annotate_scripts/orthologs_mouse2human.csv",
                input_species="mouse"):
    """transfer mouse gene names to human ortholog gene names

    Parameters
    ---------
    ref: dataframe / matrix
        mouse - human ortholog genes corresponding table
    input_names: string
        input gene names
    input_species: string
        transfer mouse to human ortholog genes or rather

    Returns
    ---------
    show and save picture

    """
    utils.warn(" --- Load Data... ---")
    genes = input_names.split(",")
    # input_data = pd.read_csv(input_file, sep=',')
    ref_file = os.path.abspath(ref)
    ref_id = {}
    ref_syb = {}
    utils.warn(" --- Create Reference Orthologues between Human and Mouse ---")
    if input_species == "mouse":
        for line in open(ref_file, "r"):
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
        for line in open(ref_file, "r"):
            # gene: human_id human_name mouse_id mouse_n
            line = line.split(",")
            h_gene_id = line[0]
            h_gene_syb = line[1]
            m_gene_id = line[2]
            m_gene_syb = line[3]
            ref_id[h_gene_id] = m_gene_id
            ref_syb[h_gene_syb] = m_gene_syb
    utils.warn(" --- Begin Transferring ---")
    for i in range(1, len(genes)):
        if genes[i] in ref_id.keys() and ref_id[genes[i]] != '':
            genes[i] = ref_id[genes[i]]
        elif genes[i] in ref_syb.keys() and ref_syb[genes[i]] != '':
            genes[i] = ref_syb[genes[i]]
    genes[0] = "Gene"

    return ','.join(genes).replace('\n', '')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Rename mouse gene name to human orthologs gene name')
    input_str = sys.stdin.read()
    parser.add_argument('--ref', '-r', required=False, help='biomart sub-database',
                        default="/home/user/data2/rbase/spatial_annotate_scripts/orthologs_mouse2human.csv")
    parser.add_argument('--input_species', '-s', required=True, help='input gene name is mouse or human')
    args = parser.parse_args()
    print(rename_gene(input_names=input_str, **vars(args)))
