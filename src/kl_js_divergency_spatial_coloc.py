import sys
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
import scipy.stats
import random
from itertools import *
import networkx as nx
from coloc_cluster_draw import *
import bootstrap_coloc_parallel as bcp


def as_dummy_df(data_input, col_cell_type='cell_type'):
    """Transform cell type long dataframe to wider table with original cell index
    Parameters:
    ------------
    data_input: dataframe
        Includes default column "cell_type" and index "cell id"
    col_cell_type: string
        Column name of cell type,  cell type names must be syntactically valid

    Returns:
    data_out: dataframe
        dummy data frame with cell types on columns
    """
    data_input["id_add"] = data_input.index
    data_input["col_add"] = data_input[col_cell_type]
    data_input["value"] = 1
    data_out = data_input[["id_add", "col_add", "value"]].pivot_table(index="id_add", columns="col_add", values="value",
                                                                      fill_value=0)
    if (data_out.index == data_input.index).all:
        return data_out
    else:
        print("incorrect cell index ")
        sys.exit(1)


def sp_grid_kern_bin(data, coord, min_num=100, h=20, n=100j, tot_num=True):
    """ construct the map of grid from data and return Kenerl density on grid

    Parameters
    ---------
    data: dataframe
        Cell type dummy table
    coord: dataframe
        Coordinates data with X and Y columns
    min_num: int
        Minimum cells count for each cell type
    h: int
        Bandwidths for x and y directions, more details in KernelDensity function in sklearn.neighbors package
    n: int + j
        Number of grid points in each direction
    tot_num: bool to decide whether to normalize

    Returns:
    ---------
    data_out: dataframe with X1,X2 coords and kernel density
        data_out with cell types columns representing kernel density in each spot.

    """

    # for cells with less than min_num, randomly add 1
    data_ = data
    if min_num < 0: min_num = 1
    col_min = data.sum() < min_num
    col_min_index = col_min[col_min].index
    if len(col_min_index) != 0:
        for i in range(0, len(col_min_index)):
            random.seed(i)
            cell_type = col_min_index[i]
            ct_dummy_array = data_.loc[:, cell_type]
            random_cell = random.sample(population=list(np.where(ct_dummy_array == 0)[0]),
                                        k=min_num - ct_dummy_array.sum())
            data_.loc[:, cell_type][random_cell] = 1
            print(col_min_index[i] + " cell randomly add to " + str(data_.loc[:, col_min_index[i]].sum()) + " cells")

    # coordinates for X1 and X2
    coord = pd.DataFrame(coord)
    coord.columns = ['X' + str(i) for i in range(0, len(coord.columns))]
    coord.index = list(data_.index)
    data_merge = pd.concat([coord, data_], axis=1)

    kde2d = KernelDensity(bandwidth=h, kernel="gaussian")
    # build 2D grid as sample points
    xx, yy = np.mgrid[data_merge['X0'].min():data_merge['X0'].max():n,
             data_merge['X1'].min():data_merge['X1'].max():n]
    xy_sample = pd.DataFrame(np.vstack([xx.ravel(), yy.ravel()]).T, columns=["X", "Y"])
    data_out = pd.DataFrame(xy_sample)
    data_out.columns = ["X0", "X1"]

    # print("estimate gaussian kernel 2D density for each cell types...")
    for i in data_.columns:
        xy_train = data_merge[["X0", "X1"]][data_merge[i] == 1]
        kde2d.fit(xy_train)
        # score_samples() returns the log-likelihood of the samples
        z = np.exp(kde2d.score_samples(xy_sample))

        # plt.figure(figsize=(6,8))
        # plt.pcolormesh(xx, yy, np.reshape(z, xx.shape), cmap="YlGnBu")
        # plt.scatter(xy_train['X0'], xy_train['X1'], s=2, facecolor="white")
        # plt.savefig("Cell type gaussian kernel density/" +i+".pdf")
        z = pd.DataFrame(z, columns=[i])
        data_out = pd.concat([data_out, z], axis=1)

    # data_out with cell types columns representing kernel density in each spot. (normalization)
    if tot_num:
        data_out = data_out.drop(["X0", "X1"], axis=1)
        data_out = data_out / data_out.sum()

    return data_out


def KL_JS_Divergence(X, eps=1e-20, diver="KL"):
    """ calculate Kullback-Leibler or Jensen-Shannon diversity

    Parameters
    ---------
    X: dataframe
        Density matrix
    eps: double flout
        small value added
    diver: str
        "KL" or "JS" representing Kullback-Leibler or Jensen-Shannon Divergence

    Returns
    ---------
    KL_D: dataframe
        KL-divergence matrix

    """

    X_ = X.fillna(eps)
    X_[X_ < eps] = eps
    X_ = X_ / X_.sum()
    n_type = len(X_.columns)
    diver_matrix = pd.DataFrame(np.zeros((n_type, n_type)))
    diver_matrix.index = X_.columns
    diver_matrix.columns = X_.columns
    # print("calculate cell types pairs " + diver + " divergence...")
    if diver == "KL":
        for i in combinations(X_.columns, 2):
            KL = scipy.stats.entropy(X_[i[0]], X_[i[1]])
            diver_matrix.loc[i[0], i[1]] = KL
            diver_matrix.loc[i[1], i[0]] = KL
    else:
        for i in combinations(X_.columns, 2):
            M = (X_[i[0]] + X_[i[1]]) / 2
            JS = 0.5 * scipy.stats.entropy(X_[i[0]], M, base=2) + 0.5 * scipy.stats.entropy(X_[i[1]], M, base=2)
            diver_matrix.loc[i[0], i[1]] = JS
            diver_matrix.loc[i[1], i[0]] = JS

    return diver_matrix


def KL_JS_boot_mst_single(dummy_df, coord_df, min_num=15, boot_n=10, prop=0.8, h=20, tot_num=True, diver="JS"):
    """  calculate KL or JS divergence and use MST to generate a tree structure by Bootstrap
         to obtain consensus cell types colocalization dissimilarity matrix

    Parameters
    ---------
    dummy_df: dataframe
        dummy data frame with cell types on columns
    coord_df: dataframe
        Coordinates data with X and Y columns
    min_num: int
        Minimum cells count for each cell type
    boot_n: int
        Number of bootstraping iteration
    prop: 0-1
        Subsample preportion
    diver: String
        use KL or JS divergence

    other see in sp_grid_kern_bin function

    Returns
    ---------
    dis_cons: dataframe
        Bootstrap KL/JS divergence
    mst_cons: dataframe
        Bootstrap MST matrix
    dis_boot_array: array
        Each Bootstrap result

    """
    coord_df = pd.DataFrame(coord_df)
    n_smp = len(dummy_df)
    n_type = len(dummy_df.columns)
    dis_boot_array = np.zeros((n_type, n_type, boot_n))
    mst_cons = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns), index=list(dummy_df.columns))
    dis_cons = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns), index=list(dummy_df.columns))
    for i in range(boot_n):
        print('---- Begin Bootstrap ' + str(i + 1) + " time ----")
        random.seed(i)
        # Bootstrap
        idx = random.sample(range(n_smp), round(n_smp * prop))
        data_boot = dummy_df.iloc[idx, :]
        coord_boot = coord_df.iloc[idx, :]
        k2d_boot = sp_grid_kern_bin(data=data_boot, coord=coord_boot, min_num=min_num, h=h, tot_num=tot_num)
        dis_boot = KL_JS_Divergence(k2d_boot, eps=1e-20, diver=diver)

        # create a graph from the adjacency matrix
        graph_boot = nx.from_pandas_adjacency(dis_boot)
        # MST
        graph_mst_boot = nx.minimum_spanning_tree(graph_boot)
        mst_boot = nx.to_pandas_adjacency(graph_mst_boot)

        # cumulate bootstrap
        dis_boot_array[:, :, i] = dis_boot
        mst_cons = mst_cons + mst_boot / boot_n
        dis_cons = dis_cons + dis_boot / boot_n

    return dis_cons, mst_cons, dis_boot_array


def network_microenv(df_adjacency, out_path, to_cpdb=True, cutoff=0.5):
    """ obtain cell types microenvironment

        Parameters
        ---------
        df_adjacency: dataframe / matrix
            divergence matrix
        out_path: string
            picture name
        cutoff: 0-1
            filter out divergence lower than cutoff percentile
        to_cpdb: bool
            whether output as CellphoneDB microenvironment file

        Returns
        ---------
        output microenvironment file

        """
    # init
    microenv = pd.Series(df_adjacency.columns, index=["Microenv_" + str(cell).replace(' ', '_')
                                                      for cell in df_adjacency.columns])
    # assign top "cutoff" percent divergence value to cutoff value
    arr_adjacency = np.array(df_adjacency).ravel()
    arr_adjacency_nonzero = []
    for i in arr_adjacency:
        if i != 0:
            arr_adjacency_nonzero = np.append(arr_adjacency_nonzero, [i])
    cutoff_v = np.percentile(np.sort(arr_adjacency_nonzero), cutoff*100)

    # find microenvironment from consensus mst
    for cell in df_adjacency.columns:
        index = "Microenv_" + str(cell).replace(' ', '_')
        # find non-zero element and correspondent cell type
        non_zero_index = df_adjacency[cell].loc[
            (df_adjacency[cell] != 0) & (df_adjacency[cell] < cutoff_v)
            ].index.values
        # add interacting cell type
        if len(non_zero_index) != 0:
            microenv[index] = np.append(cell, non_zero_index)
        else:
            microenv.drop(index, inplace=True)

    if to_cpdb:
        out_csv_df = pd.DataFrame(columns=['cell_type', 'microenvironment'])
        for k, v in microenv.items():
            v = v.reshape(len(v), 1)
            k = np.array(str(k)).repeat(len(v)).reshape(len(v), 1)
            out_csv_df = pd.concat(
                [out_csv_df, pd.DataFrame(np.hstack([v, k]), columns=['cell_type', 'microenvironment'])],
                axis=0, ignore_index=True)
        out_csv_df.to_csv(path_or_buf=out_path + "microenvironment.csv", header=True, index=False, sep=',')
        return out_csv_df
    else:
        return microenv


def spatial_cell_types_coloc(sp_data_inp, cutoff, col_cell_type="cell_type", h=20, boot_n=20, out_path="./"):
    cell_type_dummy_df = as_dummy_df(sp_data_inp.obs, col_cell_type=col_cell_type)
    dis_boot_array, dis_cons, mst_cons = bcp.KL_JS_boot_mst(dummy_df=cell_type_dummy_df,
                                                            coord_df=sp_data_inp.obsm["spatial"], h=h, boot_n=boot_n)
    network_microenv(mst_cons, cutoff=cutoff, out_path=out_path + "table/")
    divergence_clustermap(dis_cons, name="cell_types_JSD", out_path=out_path)
    network_draw(mst_cons, name="cell_types_mst_network", out_path=out_path)


if __name__ == "__main__":
    spatial_adata_annotated = sc.read('data/STW-M-Brain-Stereo-seq-1/coronal_1.bin50.adata_sp_ann.clusters.h5ad')
    spatial_cell_types_coloc(sp_data_inp=spatial_adata_annotated, col_cell_type="cell_type",
                             h=20, out_path="data/STW-M-Brain-Stereo-seq-1/out/", cutoff=0.5)
