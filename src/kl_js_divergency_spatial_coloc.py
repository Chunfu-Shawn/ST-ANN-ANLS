import matplotlib.pyplot as plt
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
import scipy.stats
import random
from itertools import *
import seaborn as sns
import networkx as nx


def as_dummy_df(data_input, col_cell_type='tangram_cell_type'):
    """Transform cell type long dataframe to wider table with original cell index
    Parameters:
    ------------
    data_input: dataframe
        Includes default column "tangram_cell_type" and index "cell id"
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
    data_merge = pd.concat([data_, coord], axis=1)

    kde2d = KernelDensity(bandwidth=h, kernel="gaussian")
    # build 2D grid as sample points
    xx, yy = np.mgrid[data_merge['X0'].min():data_merge['X0'].max():n,
             data_merge['X1'].min():data_merge['X1'].max():n]
    xy_sample = pd.DataFrame(np.vstack([xx.ravel(), yy.ravel()]).T, columns=["X", "Y"])
    data_out = pd.DataFrame(xy_sample)
    data_out.columns = ["X0", "X1"]

    print("estimate gaussian kernel 2D density for each cell types...")
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
    print("calculate cell types pairs " + diver + " divergence...")
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


def KL_JS_boot_mst(dummy_df, coord_df, min_num=15, boot_n=10, prop=0.8, h=20, tot_num=True, diver="JS"):
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
        print('---- Bootstrap ' + str(i + 1) + " time ----")
        random.seed(i)
        ## Bootstrap ##
        idx = random.sample(range(n_smp), round(n_smp * prop))
        data_boot = dummy_df.iloc[idx, :]
        coord_boot = coord_df.iloc[idx, :]
        k2d_boot = sp_grid_kern_bin(data=data_boot, coord=coord_boot, min_num=min_num, h=h, tot_num=tot_num)
        dis_boot = KL_JS_Divergence(k2d_boot, eps=1e-20, diver=diver)
        dis_boot_array[:, :, i] = dis_boot

        # create a graph from the adjacency matrix
        graph_boot = nx.from_pandas_adjacency(dis_boot)
        # MST
        graph_mst_boot = nx.minimum_spanning_tree(graph_boot)
        mst_boot = nx.to_pandas_adjacency(graph_mst_boot)

        # cumulate bootstrap
        mst_cons = mst_cons + mst_boot / boot_n
        dis_cons = dis_cons + dis_boot / boot_n

    print(dis_cons, mst_cons)
    return dis_cons, mst_cons, dis_boot_array


def divergence_heatmap(matrix, name="default_divergence", out_path='./'):
    """ plot divergence clusterd heatmap

    Parameters
    ---------
    matrix: dataframe / matrix
        divergence matrix
    name: string
        picture name

    Returns
    ---------
    show and save picture;S

    """
    matrix_log = -np.log(matrix)
    # transfre -log0 (infinite) to the max value of column
    matrix_log = matrix_log.replace(np.inf, np.nan)
    matrix_log = matrix_log.replace(np.nan, matrix_log.max())
    sns_plot = sns.clustermap(
        matrix_log,
        cmap="YlOrRd",
        linewidths=1,
        linecolor="white",
        square=True,
        cbar_pos=[.8, .55, .02, .2],
        method="ward")
    # mask upper triangle
    mask = np.triu(np.ones_like(matrix_log))
    values = sns_plot.ax_heatmap.collections[0].get_array().reshape(matrix_log.shape)
    new_values = np.ma.array(values, mask=mask)
    sns_plot.ax_heatmap.collections[0].set_array(new_values)
    # set left dendrogram invisible
    sns_plot.ax_row_dendrogram.set_visible(False)
    # set y axis ticks left
    sns_plot.ax_heatmap.yaxis.set_ticks_position("left")
    plt.show()
    # save figure
    sns_plot.savefig(out_path + name + ".pdf")


def draw_graph(df_adjacency, name="graph", node_size=20, edge_width=1, out_path='./'):
    """ plot network graph

    Parameters
    ---------
    df_adjacency: dataframe / matrix
        adjacency matrix
    name: string
        picture name
    node_size: control node size
    edge_width: control edge width
    out_path: output pathway

    Returns
    ---------
    show and save picture;

    """
    # create
    G = nx.from_pandas_adjacency(df_adjacency)
    # define position of nodes and labels
    pos = nx.planar_layout(G)
    labels_pos = {}
    for key, value in pos.items():
        random.seed(value[0])
        labels_pos[key] = (value[0], value[1])

    plt.figure(figsize=(6, 6))
    plt.axis('off')
    # plot nodes and edges of network graph
    nx.draw_networkx_nodes(G, pos=pos,
                           node_size=[1 * node_size * (item[1] + 1) for item in G.degree()],
                           label=True,
                           cmap=plt.cm.Blues)
    nx.draw_networkx_edges(G, pos=pos,
                           edge_color=[np.log(1 / d["weight"]) for (u, v, d) in G.edges(data=True)],
                           width=[np.log(1 * edge_width / d["weight"]) for (u, v, d) in G.edges(data=True)],
                           edge_cmap=plt.cm.Blues)
    nx.draw_networkx_labels(G, pos=labels_pos, font_size=6, font_weight="bold")
    # plt.show(block=False)

    # save figure
    plt.savefig(out_path + name + '.pdf', pad_inches=0.1)
    plt.show()


def spatial_cell_types_coloc(sp_data_inp, col_cell_type="tangram_cell_type", h=20, boot_n=20, out_path="./"):
    cell_type_dummy_df = as_dummy_df(sp_data_inp.obs, col_cell_type=col_cell_type)
    dis_cons, mst_cons, dis_boot_array = KL_JS_boot_mst(dummy_df=cell_type_dummy_df,
                                                        coord_df=sp_data_inp.obsm["spatial"], h=h, boot_n=boot_n)
    divergence_heatmap(dis_cons, name="cell_types_JSD", out_path=out_path)
    draw_graph(mst_cons, name="cell_types_mst_network", out_path=out_path)


if __name__ == "__main__":
    spatial_adata_annotated = sc.read('data/spatial_adata_annotated.h5ad')
    spatial_cell_types_coloc(sp_data_inp=spatial_adata_annotated, col_cell_type="tangram_cell_type",
                     h=20, boot_n=5, out_path="coloc_figures/")

