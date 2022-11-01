import torch
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
# runTangram9
import tangram as tg


def runTangram(sc_data, sp_data, marker=None, mode="cells", out="./result", key_deg="rank_genes_groups",
               use_raw_sc=None, use_raw_sp=None, cluster_label=None, top_n_markers=100):
    # Read in ad_sc & ad_sp
    warn("Read in ad_sc and ad_sp")
    ad_sc = sc_data
    ad_sp = sp_data

    if use_raw_sc:
        ad_sc = ad_sc.raw.to_adata().copy()
    if use_raw_sp:
        ad_sp = ad_sp.raw.to_adata().copy()

    # Get marker
    warn("Get marker")
    if marker is not None:
        marker = pd.read_csv(marker, header=None)
        marker = list(marker[0])
    else:
        marker = pd.DataFrame(ad_sc.uns[key_deg]['names']).head(top_n_markers)
        marker = np.array(marker).flatten().tolist()

    # Run pp_adatas
    warn("Run pp_adatas")
    tg.pp_adatas(adata_sc=ad_sc, adata_sp=ad_sp, genes=marker)
    if not ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']:
        print("Training genes in ad_sc and ad_sp are not identical!")
        sys.exit(1)

    warn("Run tangram mapping")
    # Mapping

    with torch.cuda.device(1):
        if mode == 'cells':
            ad_map = tg.map_cells_to_space(
                adata_sc=ad_sc,
                adata_sp=ad_sp,
                device="cuda",
                mode=mode
            )
        else:
            ad_map = tg.map_cells_to_space(
                adata_sc=ad_sc,
                adata_sp=ad_sp,
                device="cuda",
                mode=mode,
                cluster_label=cluster_label
            )
    torch.cuda.empty_cache()

    # Save
    warn("Save mapping result")
    ad_map.write(out + '.' + mode + '.h5ad')
    return out + '.' + mode + '.h5ad'


def warn(message):
    sys.stderr.write(message + "\n")


if __name__ == "__main__":
    descrition = 'Perform tangram mapping.'
    # 添加GPU环境变量
    os.environ['CUDA_VISIBLE_DEVICES'] = '0,1'
    os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
    adata = sc.read('data/adata.addDEG.h5ad')
    spatial_adata = sc.read('data/adata_a2p2.telen.m500.log1p.leiden.deg.h5ad')
    result_path = runTangram(adata, spatial_adata, mode="clusters", cluster_label="cell type",
                             out='data/tg_ann_result')
    res_adata = sc.read(result_path)
    tg.project_cell_annotations(res_adata, spatial_adata, annotation="cell type")
    spatial_adata.obs["tangram_cell_type"] = spatial_adata.obsm["tangram_ct_pred"].idxmax(axis=1)
    spatial_adata.write('data/spatial_adata_annotated.h5ad')
