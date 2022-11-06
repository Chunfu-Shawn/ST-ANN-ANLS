import pandas as pd
import numpy as np
import random
import networkx as nx
from multiprocessing import Pool, Manager
from kl_js_divergency_spatial_coloc import sp_grid_kern_bin, KL_JS_Divergence
import utils


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
    utils.info('### Begin Predicting Spatially Co-localized Cell Types by Bootstrap ###')
    utils.info('### Methods: ' + diver + ', MST ###')
    n_smp = len(dummy_df)
    n_type = len(dummy_df.columns)
    # Process Pool for maximum 3 processes
    pool = Pool(processes=10)
    # shared data
    manager = Manager()
    res_dict = manager.dict()
    res_dict["dis_boot_array"] = np.zeros((n_type, n_type, boot_n))
    # dis_cons
    res_dict["dis_cons"] = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns),
                                        index=list(dummy_df.columns))
    # mst_cons
    res_dict["mst_cons"] = pd.DataFrame(np.zeros((n_type, n_type)), columns=list(dummy_df.columns),
                                        index=list(dummy_df.columns))
    coord_df = pd.DataFrame(coord_df)
    for i in range(boot_n):
        random.seed(i)
        # Bootstrap
        idx = random.sample(range(n_smp), round(n_smp * prop))
        data_boot = dummy_df.iloc[idx, :]
        coord_boot = coord_df.iloc[idx, :]
        pool.apply_async(boot, (res_dict, i, boot_n, data_boot, coord_boot, min_num, h, tot_num, diver))
    # 扔了 1000个进程进进程池后，关闭进程池，不允许新的进程加入
    pool.close()
    # 运行进程池中的进程
    pool.join()

    return res_dict["dis_boot_array"], res_dict["dis_cons"], res_dict["mst_cons"]


def boot(res_dict, i, boot_n, data_boot, coord_boot, min_num, h, tot_num, diver):
    print('---- Bootstrap ' + str(i + 1) + " time Begin ----")
    k2d_boot = sp_grid_kern_bin(data=data_boot, coord=coord_boot, min_num=min_num, h=h, tot_num=tot_num)
    dis_boot = KL_JS_Divergence(k2d_boot, eps=1e-20, diver=diver)

    # create a graph from the adjacency matrix
    graph_boot = nx.from_pandas_adjacency(dis_boot)
    # MST
    graph_mst_boot = nx.minimum_spanning_tree(graph_boot)
    mst_boot = nx.to_pandas_adjacency(graph_mst_boot)
    # cumulate bootstrap
    res_dict["dis_boot_array"][:, :, i] = dis_boot
    res_dict["dis_cons"] = res_dict["dis_cons"] + dis_boot / boot_n
    res_dict["mst_cons"] = res_dict["mst_cons"] + mst_boot / boot_n
    print('---- Bootstrap ' + str(i + 1) + " time End ----")
