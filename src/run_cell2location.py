"""
File Name: run_cell2location.py
Author: Up Lee
Mail: uplee@pku.edu.cn
Created Time: Wed 26 Oct 2022 11:18:34 AM CST
"""

import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel, Cell2location


def read_sp(input, use_raw=True):
    adata = sc.read_h5ad(input)

    # use raw count
    if use_raw:
        adata = adata.raw.to_adata()

    # remove mito genes
    if 'mt' in adata.var.columns:
        adata = adata[:, ~adata.var['mt']].copy()

    adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata


def read_sc(count, meta, celltype_key='cell_type', batch_key=None, categorical_covariate_key=None,
            continuous_covariate_key=None):
    count = sc.read_csv(count)
    count.X = scipy.sparse.csc_matrix(count.X)

    meta = pd.read_csv(meta, index_col=0)

    if all(count.obs.index == meta.index):
        count.obs = pd.concat([count.obs, meta], axis=1)

    if celltype_key not in count.obs.columns:
        raise KeyError('celltype_key ' + celltype_key + ' not found in ' + ', '.join(count.obs.columns) + '!')
    if batch_key != None and batch_key not in count.obs.columns:
        raise KeyError('batch_key ' + batch_key + ' not found in ' + ', '.join(count.obs.columns) + '!')
    if categorical_covariate_key != None and not all(x in count.obs.columns for x in categorical_covariate_key):
        raise KeyError('categorical_covariate_key ' + categorical_covariate_key + 'is not a subset of ' + ', '.join(
            count.obs.columns) + '!')
    if continuous_covariate_key != None and not all(x in count.obs.columns for x in continuous_covariate_key):
        raise KeyError('continuous_covariate_key ' + continuous_covariate_key + 'is not a subset of ' + ', '.join(
            count.obs.columns) + '!')

    return count


class PerformCell2loc:
    def __init__(self, adata_sp, adata_sc, output_dir, celltype_key='cell_type', batch_key='batch',
                 categorical_covariate_key=None, continuous_covariate_key=None, N_cells_per_location=30,
                 detection_alpha=20):
        #### Filter
        adata_sp = adata_sp[:, ~adata_sp.var['mt']].copy()

        selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
        adata_sc = adata_sc[:, selected].copy()

        self.adata_sp = adata_sp
        self.adata_sc = adata_sc
        self.celltype_key = celltype_key
        self.batch_key = batch_key
        self.categorical_covariate_key = categorical_covariate_key
        self.continuous_covariate_key = continuous_covariate_key
        self.N_cells_per_location = N_cells_per_location
        self.detection_alpha = detection_alpha
        self.output_dir = output_dir

        print("adata_sc self int: ", self.adata_sc)

    def NB_regression(self):
        adata_sc = self.adata_sc

        # prepare anndata for the regression model
        RegressionModel.setup_anndata(adata=adata_sc,
                                      # 10X reaction / sample / batch
                                      batch_key=self.batch_key,
                                      # cell type, covariate used for constructing signatures
                                      labels_key=self.celltype_key,
                                      # multiplicative technical effects (platform, 3' vs 5', donor effect)
                                      categorical_covariate_keys=self.categorical_covariate_key,
                                      continuous_covariate_keys=self.continuous_covariate_key
                                      )

        # create the regression model
        mod = RegressionModel(adata_sc)

        # view anndata_setup as a sanity check
        mod.view_anndata_setup()

        # train mod
        mod.train(max_epochs=250, use_gpu=True)
        mod.plot_history(20)

        # Estimated cell abundance (summary of the posterior distribution)
        adata_sc = mod.export_posterior(
            adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
        )

        # QC
        mod.plot_QC()

        # Export estimated expression in each cluster
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                             for i in adata_sc.uns['mod']['factor_names']]].copy()

        inf_aver.columns = adata_sc.uns['mod']['factor_names']

        print("cell type signature: ", inf_aver)
        print("adata_sc: ", adata_sc)
        print("adata_sc self: ", self.adata_sc)

        self.mod_sc = mod
        self.signature = inf_aver

    def spatial_mapping(self):
        adata_sc = self.adata_sc
        adata_sp = self.adata_sp
        signature = self.signature

        ## Find shared genes
        # spatial transcriptomic
        intersect = np.intersect1d(adata_sp.var_names, signature.index)
        adata_sp = adata_sp[:, intersect].copy()
        # cell type signature
        signature = signature.loc[intersect, :].copy()

        ## Construct model
        # prepare anndata for cell2location model
        Cell2location.setup_anndata(adata=adata_sp, batch_key=None)

        ## Create and train the model
        mod = Cell2location(
            adata_sp, cell_state_df=signature,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=self.N_cells_per_location,
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            detection_alpha=self.detection_alpha
        )

        mod.view_anndata_setup()
        mod.train(max_epochs=10000,
                  # train using full data (batch_size=None)
                  batch_size=None,
                  # use all data points in training because
                  # we need to estimate cell abundance at all locations
                  train_size=1,
                  use_gpu=True)
        mod.plot_history(1000)

        ## Estimated cell abundance (summary of the posterior distribution)
        adata_sp = mod.export_posterior(
            adata_sp, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
        )

        self.adata_sp = adata_sp

        ## QC
        mod.plot_QC()

        self.mod_sp = mod

    def recapitulate_abundance(self):
        adata_sp = self.adata_sp
        print("adata_sp: ", adata_sp)

        adata_sp.obs[adata_sp.uns['mod']['factor_names']] = adata_sp.obsm['q05_cell_abundance_w_sf']

    def save_results(self):
        adata_sc = self.adata_sc
        adata_sp = self.adata_sp
        mod_sc = self.mod_sc
        mod_sp = self.mod_sp
        output_dir = self.output_dir

        # Save NB regression model
        mod_sc.save(output_dir, overwrite=True)
        os.rename(os.path.join(output_dir, 'model.pt'), os.path.join(output_dir, 'model.sc.pt'))

        # Save spatial mapping model
        mod_sp.save(os.path.join(output_dir), overwrite=True)
        os.rename(os.path.join(output_dir, 'model.pt'), os.path.join(output_dir, 'model.sp.pt'))

        # Save scRNA-seq anndata with cell type signature
        adata_sc.write_h5ad(os.path.join(output_dir, 'adata_sc.cellTypeSignature.h5ad'))

        # Save spatial transcriptome anndata with spatial mapping
        adata_sp.write_h5ad(os.path.join(output_dir, 'adata_sp.spatialMapping.h5ad'))


def main(argsv):
    parser = argparse.ArgumentParser(description="Process Spatial transcriptome data.")
    parser.add_argument("--dataset", type=str, help="STW id")
    parser.add_argument("--section", type=str, help="The section id of the dataset.")
    parser.add_argument("--label", type=str, default="", help="The unique identifier in name of input and output")
    parser.add_argument('--input_path', type=str, help='The input directory of process_st sp data')
    parser.add_argument('--sc_count', type=str, help='The path to count matrix of scRNA-seq')
    parser.add_argument('--sc_meta', type=str, help='The path to meta table of scRNA-seq')
    parser.add_argument("--output_path", type=str, help='The output directory of results')
    parser.add_argument("--celltype_key", type=str, default="cell_type",
                        help='The column name indicating cell type in scRNA-seq metadata table')
    parser.add_argument("--batch_key", type=str, default=None,
                        help='The column name indicating batch effect in scRNA-seq metadata table')
    parser.add_argument("--categorical_covariate_key", type=str, nargs="+", default=None,
                        help='The column name indicating categorical covariate, this argument can be passed in a continuous way')
    parser.add_argument("--continuous_covariate_key", type=str, nargs="+", default=None,
                        help='The column name indicating continuous covariate, this argument can be passed in a continuous way')
    parser.add_argument("--N_cells_per_location", type=int, default=30, help='N_cells_per_location')
    parser.add_argument("--detection_alpha", type=int, default=20, help='detection_alpha')
    parser.add_argument("--uid", type=str, help='Unique ID')
    parser.add_argument("--gpu_index", type=int, default=0, choices=[0, 1], help='The index of GPU to use')
    args = parser.parse_args()

    # print("covariate_key: ",args.covariate_key)

    if args.label != "":
        input = os.path.join(args.input_path, args.dataset, args.section,
                             '.'.join([args.section, args.label, "process_st.h5ad"]))
    else:
        input = os.path.join(args.input_path, args.dataset, args.section, '.'.join([args.section, "process_st.h5ad"]))

    output_dir = os.path.join(args.output_path, args.dataset, args.section, '_'.join(['cell2location', args.uid]))
    os.makedirs(output_dir, exist_ok=True)

    ######################################
    #### Read data
    ######################################
    print("#### Read data")
    adata_sp = read_sp(input=input, use_raw=True)
    adata_sc = read_sc(args.sc_count,
                       args.sc_meta,
                       celltype_key=args.celltype_key,
                       batch_key=args.batch_key,
                       categorical_covariate_key=args.categorical_covariate_key,
                       continuous_covariate_key=args.continuous_covariate_key)

    print("adata_sp: ", adata_sp)
    print("adata_sc: ", adata_sc)

    ######################################
    #### Cell2location
    ######################################
    print("#### Estimate cell type signature using NB regression")
    cell2loc = PerformCell2loc(adata_sp=adata_sp,
                               adata_sc=adata_sc,
                               output_dir=output_dir,
                               celltype_key=args.celltype_key,
                               batch_key=args.batch_key,
                               categorical_covariate_key=args.categorical_covariate_key,
                               continuous_covariate_key=args.continuous_covariate_key,
                               N_cells_per_location=args.N_cells_per_location,
                               detection_alpha=args.detection_alpha
                               )
    print("#### NB_regression")
    cell2loc.NB_regression()
    print("#### Spatial mapping")
    cell2loc.spatial_mapping()
    print("#### Retrieve cell abundance")
    cell2loc.recapitulate_abundance()
    print("#### Save results")
    cell2loc.save_results()


if __name__ == '__main__':
    import sys

    main(sys.argv)
