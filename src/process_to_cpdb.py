import numpy as np
import pandas as pd
import scanpy as sc
import gene_orthologues_trans as got
import utils


def process_to_cpdb(sc_data, save_path, name="test", annotation_name="cell_type"):
    adata_sc.X = sc_data.X.A

    utils.warn(' --- Generating count file --- ')
    # set df.dtypes = float32
    df_expr_matrix = pd.DataFrame(sc.pp.normalize_total(adata_sc, inplace=False)["X"], dtype="float")
    df_expr_matrix = df_expr_matrix.T
    # set cell ids as columns
    df_expr_matrix.columns = adata_sc.obs.index
    # genes should be either Ensembl IDs or gene names
    # transform orthologs
    index = adata_sc.var.index.values
    genes = got.rename_gene(np.array(index))
    df_expr_matrix.set_index(genes, inplace=True)

    utils.warn(' --- Saving count file... --- ')
    df_expr_matrix.to_csv(save_path + name + ".counts.txt", sep='\t')

    utils.warn(' --- Generating metadata file --- ')
    df_meta = pd.DataFrame(
        data={'Cell': list(adata_sc.obs.index), 'cell_type': list(adata_sc.obs[annotation_name])}
    )
    df_meta.set_index('Cell', inplace=True)
    utils.warn(' --- Saving metadata file... --- ')
    df_meta.to_csv(save_path + name + ".meta.txt", sep='\t')


if __name__ == "__main__":
    adata_sc = utils.read_sc(
        count='../data/STW-M-Brain-Stereo-seq-1/counts.csv.gz',
        meta='../data/STW-M-Brain-Stereo-seq-1/labels.csv.gz',
        celltype_key="cell_type"
    )
    process_to_cpdb(
        sc_data=adata_sc,
        annotation_name="cell_type",
        save_path="../data/STW-M-Brain-Stereo-seq-1/",
        name="coronal_1"
    )



