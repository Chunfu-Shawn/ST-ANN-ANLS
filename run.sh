################################################
#File Name: run.sh
#Author: rbase    
#Mail: uplee@pku.edu.cn
#Created Time: Tue 28 Jun 2022 11:42:33 AM CST
################################################

#!/bin/sh 

scriptDir=/home/user/data2/uplee/projects/spatialTransWeb/bin
ad_sc=/home/user/data2/uplee/projects/spatialTransWeb/data/public/scRNAseq/mouseCortex/data/annData_annotation_qijt/adata.addDEG.h5ad
ad_sp=/home/user/data2/uplee/projects/spatialTransWeb/spatial/inhouse/apcdd1_e14.5_brain/add_gem/bin20/spatial_cluster/outs/A2-2/adata_a2p2.telen.m500.log1p.leiden.deg.h5ad
out_put=/home/user/data2/rbase/test_tangram


time (python $scriptDir/run_tangram_mapping.py \
  --ad_sc $ad_sc \
  --ad_sp $ad_sp \
  --use_raw_sc \
  --use_raw_sp \
  --key_deg rank_genes_groups_ct \
  --top_n_marker 100 \
  --device cuda \
  --mode cells \
  --cluster_label 'cell type' \
  --out $out_put/admap_clsCt_a2p2.telen.m500.log1p.leiden.deg.h5ad) 1>$out_put/tangram.a2p2_telen.clsCt.log 2>&1
