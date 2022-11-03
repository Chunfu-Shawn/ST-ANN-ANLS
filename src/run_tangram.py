import torch
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import utils
import argparse
# runTangram
import tangram as tg


class PerformTangram:
    def __init__(self, ad_sc, ad_sp, output_dir='./', cell_type_key='cell_type', use_raw_sc=None,
                 use_raw_sp=None, mode="clusters"):
        # Read in ad_sc & ad_sp
        utils.warn("Read in ad_sc and ad_sp")
        self.adata_sc = ad_sc
        self.adata_sp = ad_sp[:, ~adata_sp.var['mt']].copy()
        if use_raw_sc:
            self.adata_sc = ad_sc.raw.to_adata().copy()
        if use_raw_sp:
            self.adata_sc = ad_sp.raw.to_adata().copy()

        self.cell_type_key = cell_type_key
        self.output_dir = output_dir
        self.mode = mode
        self.adata_map = None
        self.adata_sp_ann = None
        self.adata_ge = None

    def run_tangram(self, marker=None, key_deg="rank_genes_groups", top_n_markers=100, gpu_index=0):
        # Get marker
        utils.warn("Get marker")
        if marker is not None:
            marker = pd.read_csv(marker, header=None)
            marker = list(marker[0])
        else:
            marker = pd.DataFrame(self.adata_sc.uns[key_deg]['names']).head(top_n_markers)
            marker = np.array(marker).flatten().tolist()

        # Run pp_adatas
        utils.warn("Run pp_adatas")
        tg.pp_adatas(adata_sc=self.adata_sc, adata_sp=self.adata_sp, genes=marker)
        if not self.adata_sc.uns['training_genes'] == self.adata_sp.uns['training_genes']:
            print("Training genes in ad_sc and ad_sp are not identical!")
            sys.exit(1)

        utils.warn("Run tangram mapping")
        # Mapping

        with torch.cuda.device(gpu_index):
            if self.mode == 'cells':
                adata_map = tg.map_cells_to_space(
                    adata_sc=self.adata_sc,
                    adata_sp=self.adata_sp,
                    device="cuda",
                    mode=self.mode
                )
            else:
                adata_map = tg.map_cells_to_space(
                    adata_sc=self.adata_sc,
                    adata_sp=self.adata_sp,
                    device="cuda",
                    mode=self.mode,
                    cluster_label=self.cell_type_key
                )
        torch.cuda.empty_cache()
        self.adata_map = adata_map

    def project_cell_annotations(self):
        utils.warn("Project prediction cell types to ST data")
        tg.project_cell_annotations(self.adata_map, self.adata_sp, annotation=self.cell_type_key)
        adata_sp.obs[self.cell_type_key] = self.adata_sp.obsm["tangram_ct_pred"].idxmax(axis=1)
        self.adata_sp_ann = adata_sp

    def project_genes(self):
        self.adata_ge = tg.project_genes(self.adata_map, self.adata_sc, cluster_label=self.cell_type_key)

    def return_adata_sp_ann(self):
        return self.adata_sp_ann

    def save(self):
        utils.warn("Save tangram results")
        if self.adata_map:
            self.adata_map.write(self.output_dir + 'adata_map.' + self.mode + '.h5ad')
        if self.adata_sp_ann:
            self.adata_sp_ann.write(self.output_dir + 'adata_sp_ann.' + self.mode + '.h5ad')
        if self.adata_ge:
            self.adata_ge.write(self.output_dir + 'adata_ge.' + self.mode + '.h5ad')


if __name__ == "__main__":
    utils.warn("Perform tangram mapping...")
    # 添加GPU环境变量
    os.environ['CUDA_VISIBLE_DEVICES'] = '0,1'
    os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
    # sc_adata = utils.read_sc('data/adata.addDEG.h5ad')
    # sp_adata = utils.read_sp('data/adata_a2p2.telen.m500.log1p.leiden.deg.h5ad')
    adata_sc = sc.read('data/adata.addDEG.h5ad')
    adata_sp = sc.read('data/adata_a2p2.telen.m500.log1p.leiden.deg.h5ad')
    tg_ins = PerformTangram(ad_sc=adata_sc,
                            ad_sp=adata_sp,
                            mode="clusters",
                            output_dir='data/',
                            cell_type_key='cell type')
    tg_ins.run_tangram(gpu_index=1)
    tg_ins.project_cell_annotations()
    print(tg_ins.return_adata_sp_ann())
