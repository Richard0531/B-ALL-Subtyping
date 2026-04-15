library(DESeq2)
library(Rtsne)
library(pals)
library(uwot)
library(tidyverse)
library(dplyr)
library(sva)
library(edgeR)
library(colorspace)
library(readxl)
library(MDALL)
library(ggrepel)
args <- commandArgs(trailingOnly = TRUE)
setwd('/home/hcliulab/yhsiao/gene_expression/')
set.seed(123)

sample_name = args[1]
dir_name = args[2]
#dff <- read.table(paste0(dir_name,sample_name,'_count_multi_simplified.txt'), header = TRUE, row.names = 1)
dff <- read.table(paste0(dir_name,'/',sample_name,'/',sample_name,'_count_multi_simplified.txt'), header = TRUE, row.names = 1)
colnames(dff)[1] = 'sample'
case_list = colnames(dff)[1] 

dff <- data.frame(feature = rownames(dff),dff,check.names = FALSE)
df_vst=get_vst_values(obj_in = obj_234_HTSeq,df_count = dff)
df_vst_i=f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vst)
obj_=obj_merge(obj_in = obj_1821,df_in = df_vst_i,assay_name_in = "vst")
variable_genes_in=obj_[['keyFeatures']][1:min(c(1058),length((obj_[['keyFeatures']])))]
fit_umap=umap(t(assays(obj_$SE[variable_genes_in,])[["vst"]]),n_neighbors=10,n_components=2,min_dist=0.3)
df_umap=as.data.frame(fit_umap)
names(df_umap)=c("uMAP_1","uMAP_2","uMAP_3")[1:ncol(df_umap)]
df_umap$n_neighbors=10
df_umap$FeatureN=dim(t(assays(obj_$SE[variable_genes_in,])[["vst"]]))[2]
df_umap$COH_sample=row.names(df_umap)
obj_[['umap']]=df_umap
obj_[[paste0('umap',"_fit")]]=fit_umap
df_plot=get_embeding_feature(obj_in = obj_,feature = 'diag_raw',reduction = 'umap')
p1=gg_dimPlot(df_in=df_plot,x="uMAP_1",y="uMAP_2",var_col='diag_raw',cols=subtypeCol(),
              size=1,axis_by=10,
              highlightLevel='TestSample',sizeHighlight=2,
              plot_title='UMAP',x_lab="UMAP1",y_lab="UMAP2",legend_title='Subtypes')

ref.subtypes <- df_plot$diag_raw
t1.mean <- tapply(df_plot[1:length(ref.subtypes),2], ref.subtypes, mean)
t2.mean <- tapply(df_plot[1:length(ref.subtypes),3], ref.subtypes, mean)
t2.mean = t2.mean+2.5
df_labels = data.frame(t1.mean,t2.mean,label = rownames(t1.mean))
df_labels = df_labels[df_labels$label != 'TestSample',]
df_labels[1,2] = df_labels[1,2] -1

pdf(paste0(dir_name,'/',sample_name,'/','case_gene_expression_BALL_MDALL.pdf'), width = 8, height = 8)
p1 + geom_text(
  data = df_labels,
  aes(x = t1.mean, y = t2.mean, label = label),
  inherit.aes = T,
  size = 4
)
dev.off()































