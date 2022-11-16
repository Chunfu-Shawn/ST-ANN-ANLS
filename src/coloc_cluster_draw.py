import seaborn as sns
import networkx as nx
from scipy.cluster.hierarchy import ward, fcluster
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json


def output2echarts_divergence_cluster(matrix, name, out_path="./"):
    """ plot divergence clusterd heatmap

    Parameters
    ---------
    matrix: dataframe / matrix
        divergence matrix
    out_path: string
        output directory


    Returns
    ---------
    cluster matrix

    """
    # precompute linkage matrix
    Z = ward(matrix)
    index_begin = len(matrix.index)
    rearrange = {}
    for i in range(Z.shape[0]):
        x = [Z[i, 0]] if Z[i, 0] <= Z.shape[0] else rearrange[int(Z[i, 0])]
        y = [Z[i, 1]] if Z[i, 1] <= Z.shape[0] else rearrange[int(Z[i, 1])]
        x.extend(y)
        rearrange[index_begin] = x
        index_begin += 1
    cluster_list = rearrange[list(rearrange.keys())[-1]]
    cluster_matrix = matrix.iloc[cluster_list, cluster_list]
    cluster_matrix.index = cluster_matrix.columns
    # save as csv
    cluster_matrix.to_csv(path_or_buf=out_path + name + ".csv", header=True, index=False, sep=',')
    # save as json
    cell_types = list(cluster_matrix.columns)
    log_jsd = []
    for i in range(len(cell_types)):
        for j in range(len(cell_types)):
            log_jsd.append([i, j, cluster_matrix.iloc[i, j]])
    out_dict = {
        "cell_types": cell_types,
        "log_jsd": log_jsd
    }
    file = open(out_path + name + ".json", "w")
    json.dump(out_dict, file)
    file.close()
    return cluster_matrix


def output2echarts_dict_graph(G, out_path="./", name="mst"):
    nodes = []
    edges = []
    for i in list(G.degree()):
        nodes.append({"id": i[0], "name": i[0], "degree": i[1]})
    for i in list(G.edges(data=True)):
        edges.append({"source": i[0], "target": i[1], "weight": -np.log2(i[2]["weight"])})
    out_dict = {
        "nodes": nodes,
        "edges": edges
    }
    file = open(out_path + name + ".json", "w")
    json.dump(out_dict, file)
    file.close()
    return out_dict


def divergence_clustermap(matrix, name="default_divergence", out_path='./'):
    """ plot divergence clusterd heatmap

    Parameters
    ---------
    matrix: dataframe / matrix
        divergence matrix
    name: string
        picture name
    out_path: string
        output directory

    Returns
    ---------
    show and save picture

    """
    matrix_log = -np.log2(matrix)
    # transfre -log0 (infinite) to the max value of column
    matrix_log = matrix_log.replace(np.inf, np.nan)
    matrix_log = matrix_log.replace(np.nan, matrix_log.stack().max())
    # save csv and json
    output2echarts_divergence_cluster(matrix_log, name, out_path=out_path + "table/")
    # draw
    sns_plot = sns.clustermap(
        matrix_log,
        cmap="YlOrRd",
        cbar_pos=[.8, .55, .02, .2],
        dendrogram_ratio=0.1,
        method="ward")
    # mask upper triangle
    # mask = np.triu(np.ones_like(matrix_log))
    # values = sns_plot.ax_heatmap.collections[0].get_array().reshape(matrix_log.shape)
    # new_values = np.ma.array(values, mask=mask)
    # sns_plot.ax_heatmap.collections[0].set_array(new_values)
    # set left dendrogram invisible
    sns_plot.ax_row_dendrogram.set_visible(False)
    # set y axis ticks left
    sns_plot.ax_heatmap.yaxis.set_ticks_position("left")
    plt.show()
    # save figure
    sns_plot.savefig(out_path + "pdf/" + name + ".pdf")


def network_draw(df_adjacency, name="graph", node_size=20, edge_width=1, out_path='./'):
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
    show and save picture

    """
    # create
    G = nx.from_pandas_adjacency(df_adjacency)
    # define position of nodes and labels
    pos = nx.planar_layout(G)
    labels_pos = {}
    for key, value in pos.items():
        labels_pos[key] = (value[0], value[1])

    plt.figure(figsize=(6, 6))
    plt.axis('off')
    # save as csv
    df_adjacency.to_csv(out_path + "table/" + name + '.csv')
    # save as json
    output2echarts_dict_graph(G, out_path + "table/", name)
    # plot nodes and edges of network graph
    nx.draw_networkx_nodes(G, pos=pos,
                           node_size=[1 * node_size * (item[1] + 1) for item in G.degree()],
                           label=True, node_color=plt.cm.YlOrRd(0.6))
    nx.draw_networkx_edges(G, pos=pos,
                           edge_color=[-np.log2(d["weight"]) for (u, v, d) in G.edges(data=True)],
                           width=[-np.log2(d["weight"]) * edge_width for (u, v, d) in G.edges(data=True)],
                           edge_cmap=plt.cm.YlOrRd)
    nx.draw_networkx_labels(G, pos=labels_pos, font_size=6, font_weight="bold")
    # plt.show(block=False)

    # save figure
    plt.savefig(out_path + "pdf/" + name + '.pdf', pad_inches=0.1)
    plt.show()
