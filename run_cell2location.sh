################################################
#File Name: run_cell2location.sh
#Author: Up Lee    
#Mail: uplee@pku.edu.cn
#Created Time: Wed 26 Oct 2022 11:18:34 AM CST
################################################

#!/bin/sh 

#### 2022-10-26 ####

# Path to perform integration analysis

#### Test analysis using example data

# See 

#### EO 2022-10-26 ####

#### 2022-10-30 ####

###########################################
#### Run cell2location
###########################################

mkdir -p log

##################################
#### STW-M-Brain-Stereo-seq-1
 ### coronal_1
  ## bin100
   # Mouse-corticogenesis
##################################

scriptDir=/home/user/data3/uplee/projects/spatialTransWeb/bin
dataset=STW-M-Brain-Stereo-seq-1
section=coronal_1
label=bin100
inputPath=/home/user/data3/uplee/projects/spatialTransWeb/spatial/visualization/
scCount=/home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/Mouse-corticogenesis/data/counts.csv.gz
scMeta=/home/user/data3/uplee/projects/spatialTransWeb/scRNA-seq/preparation/Mouse-corticogenesis/data/labels.csv.gz
outputPath=/home/user/data3/uplee/projects/spatialTransWeb/spatial/integration
batchKey=batch
celltypeKey=cell_type

uid=$( python $scriptDir/generate_uid.py )
python $scriptDir/run_cell2location.py \
 --dataset $dataset \
 --section $section  \
 --label $label  \
 --input_path $inputPath \
 --sc_count $scCount \
 --sc_meta $scMeta \
 --output_path $outputPath \
 --batch_key $batchKey \
 --celltype_key $celltypeKey \
 --uid $uid 1>log/run_cell2location.$dataset.$section.$label.$uid.log 2>log/run_cell2location.$dataset.$section.$label.$uid.err
 
#### EO 2022-10-30 ####

