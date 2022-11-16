# author: "Chunfu-Shawn"
# date: "2022-11-04"

  
# install essential packages
library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(latex2exp)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(optparse)
library(rjson)
suppressPackageStartupMessages(library(jsonlite))


# select cell pairs from microenvironment
select_cell_pairs_from_microenv = function(microenvs_path){
  print(" --- select cell pairs from microenvironment --- ")
  microenvs_data = read_csv(microenvs_path,col_names = T, show_col_types = F)
  microenvs = unique(microenvs_data[["microenvironment"]])
  microenvs_cell_pairs = list()
  for (env in microenvs){
    cell_types = microenvs_data[microenvs_data['microenvironment'] == env,][["cell_type"]]
    for (i in cell_types){
      for (j in cell_types){
        microenvs_cell_pairs[[env]] = c(microenvs_cell_pairs[[env]],paste(i,j,sep = '|'))
      }
    }
  }
  return(microenvs_cell_pairs)
}


# ggplot dotplot drawing
dotplot = function(data_ip,output_path,output_name){
  #adjusted parameters
  xlen = length(unique(data_ip$interacting_pair))
  ylen = length(unique(data_ip$cell_pairs))
  width = xlen/15+10
  xsize = -0.02*xlen + 9
  height = ylen/12+6
  ysize = -0.08*ylen + 10
  colorlimit = quantile(data_ip$means,probs = c(0,1))
  # draw
  ggplot(data=data_ip, aes(x=interacting_pair,y=cell_pairs)) +
    geom_point(aes(size=pvalue,colour=means)) +
    scale_size_continuous(name = TeX('$-log_{10}$ p-value'), range = c(0.1,2))+ 
    scale_colour_gradientn(colours = viridis(100,option = "D"),
                           values = c(seq(0,0.3,length.out = 70),seq(0.3,1,length.out = 30)),
                           name = TeX('$log_{2}$ mean expr (molecule 1,molecule 2)')) +
    xlab('')+ylab('')+
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,size=xsize),
          axis.text.y = element_text(angle = 30,hjust = 1,size=ysize),
          legend.key = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  ggsave(paste(output_path,'/pdf/', output_name,'.pdf',sep = ''), device = 'pdf', width = width, height = height)
}


# define function to draw interaction counts heatmap 
cells_interaction_heatmap = function(
    data_ip,microenvs_path,output_path,output_name="inter_count_heatmap"){
  
  microenvs_data = read_csv(microenvs_path,col_names = T,show_col_types = FALSE)
  cell_types = unique(microenvs_data[["cell_type"]])
  count_mat = matrix(data = 0,nrow = length(cell_types),ncol = length(cell_types))
  rownames(count_mat) = cell_types
  colnames(count_mat) = cell_types
  for (i in 1:length(rownames(data_ip))){
    cell_type = strsplit(data_ip[["cell_pairs"]][i], '\\|')[[1]]
    count_mat[cell_type[1],cell_type[2]] = count_mat[cell_type[1],cell_type[2]] +1
    count_mat[cell_type[2],cell_type[1]] = count_mat[cell_type[2],cell_type[1]] +1
  }
  
  # rescale color
  count_array = c(count_mat[lower.tri(count_mat)])
  col_func = circlize::colorRamp2(quantile(count_array, probs = c(seq(0,0.7,length.out=30),seq(0.7,1,length.out=70))), viridis(100,option = "D"))
  # draw
  pdf(file = paste(output_path,'/pdf/',output_name,'.pdf',sep = ''))
  hm = Heatmap(count_mat, name = 'Number of interaction',
               col = col_func,
               heatmap_legend_param = list(
                 direction = "horizontal"
               ),
               row_names_gp = gpar(fontsize=8),
               column_names_gp = gpar(fontsize=8),
               width = unit(10,"cm"), height = unit(10,"cm"))
  draw(hm , heatmap_legend_side="bottom")
  dev.off()
  
  # cluster
  hc = hclust(dist(count_mat),"ward.D")
  
  # save as json
  count_mat = count_mat[hc$order,hc$order]
  rownames(count_mat) = NULL
  colnames(count_mat) = NULL
  count_df = expand_grid(x=seq_len(nrow(count_mat)),y=seq_len(ncol(count_mat))) %>%
    rowwise() %>%
    mutate(count=count_mat[x,y])
  count = apply(count_df, 1, function(x){
    c(as.integer(x["x"])-1,as.integer(x["y"])-1,as.integer(x["count"]))
  },simplify = F)
  out_json =toJSON(list(cell_types=hc$labels[hc$order],count=count),pretty=F,auto_unbox=T)
  cat(out_json, file = paste(output_path,'/table/',output_name,'.json',sep = ''))
}


# data filter and plot
interaction_plot = 
  function(means_path, pvalues_path,
           output_path="../data/STW-M-Brain-Stereo-seq-1/out/",
           output_name="dot_plot", microenvs_path=NULL){
    means = read_tsv(means_path, col_names = T, show_col_types = FALSE)
    pvalues = read_tsv(pvalues_path, col_names = T, show_col_types = FALSE)
    print(" --- interacting ligands and receptors data filtering... --- ")
    # cluster by interacting pair
    clust = hclust(dist(pvalues[,-(1:11)] %>% as.matrix()))
    
    # wide table to long table and log2/10
    cell_pairs = colnames(means)
    means_all = pivot_longer(means,cols = colnames(means)[-(1:11)],names_to = "cell_pairs", values_to = "means") %>%
      select(interacting_pair,cell_pairs,means) %>% 
      mutate(interacting_pair = factor(interacting_pair, levels = pvalues[clust$order,]$interacting_pair))
    means_all$means = log2(means_all$means)
    means_all$means[is.infinite(means_all$means)|means_all$means< -2] = -2
    
    pvalues_all = pivot_longer(pvalues,cols = colnames(pvalues)[-(1:11)],names_to = "cell_pairs", values_to = "pvalue") %>%
      select(interacting_pair,cell_pairs,pvalue) %>% 
      mutate(interacting_pair = factor(interacting_pair, levels = pvalues[clust$order,]$interacting_pair))
    pvalues_all[["pvalue"]] = -log10(pvalues_all[["pvalue"]])
    pvalues_all$pvalue[is.infinite(pvalues_all$pvalue)] = 6
    
    
    # filter out pvalue > 0.01 and merge pvalues and means ( follow pvalue )
    pvalues_means_all = pvalues_all %>% filter(pvalue>= 2) %>% 
      left_join(means_all, by = c("cell_pairs","interacting_pair"))
    # dot plot
    print(" --- Drawing p-value and mean expression data dotplot... --- ")
    pvalues_means_all %>% dotplot(output_path=output_path,output_name=output_name)
    
    if (!is.null(microenvs_path)){
      # interaction count heatmap
      print(" --- Drawing Number of interaction heatmap... --- ")
      cells_interaction_heatmap(data_ip = pvalues_means_all, 
                                microenvs_path=microenvs_path,
                                output_path=output_path)
      microenvs_cell_pairs = select_cell_pairs_from_microenv(
        microenvs_path=microenvs_path
      )
      ## for each environment
      print(" --- Drawing dotplot for each environment... --- ")
      for (env in names(microenvs_cell_pairs)){
        pvalues_means_env = pvalues_means_all[
          pvalues_means_all$cell_pairs %in% microenvs_cell_pairs[[env]],]
        pvalues_means_env %>% 
          dotplot(output_path=output_path,output_name=paste(output_name,env,sep = "_"))
      }
    }
  }


# operate shell parse
option_list = list(
  make_option(c("--work_dir","-w"), type = "character", default = ".", help="path of input files and output graphs")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
work_dir = opt$work_dir
work_dir = gsub("/$","",work_dir)

# set paramters
means_path = paste(work_dir,"/means.txt",sep = "")
pvalues_path = paste(work_dir,"/pvalues.txt",sep = "")
microenvs_path = paste(work_dir,"/table/microenvironment.csv",sep = "")

# main 
interaction_plot(
  means_path = means_path,
  pvalues_path = pvalues_path,
  microenvs_path = microenvs_path,
  output_path = work_dir,
)