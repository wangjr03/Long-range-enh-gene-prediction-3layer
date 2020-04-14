args <- commandArgs(T)

ct <- as.numeric(args[1])

load(paste0("../output_exons_updated/",ct,"/info_frame.Rdata"))

gene_annotation <- read.table("../gene_annotation/gene_annotation_V19.txt")

promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")

enh_coord <- read.table("../consensus_enhancer_coord.txt")


# pool <- subset(info_frame,info_frame$corrected_p>0.85&info_frame$corr>0.3&(info_frame$Thymus_CaptureC==1|info_frame$Hi_C==1)&info_frame$tf_count>0&(info_frame$GTEX!=0|info_frame$muther!=0|info_frame$battle!=0))
# pool$gene_name <- gene_annotation[match(promoter_frame[pool[,1],1],gene_annotation[,7]),6]


#TF should also be meaningful TF: enrich in gene group

select_example <- function(pool){
  
  #get gene group TF profile
  enh_coord <- read.table("../consensus_enhancer_coord.txt")
  
  
  load(paste0("../output_exons_updated/",ct,"/i30_D5_DD10_C1_G200_K0.07_T10_B1_qt_pearson_NA/iteration_55.Rdata"))
  
  
  #normalize group_TF_bar with z score
  group_TF_bar <- sapply(group_tf_data,function(x){apply(x,2,mean)})
  
  group_TF_bar <- apply(group_TF_bar,2,function(x){(x-mean(x))/sd(x)})
  
  
  #load enh-TF matrix
  
  load(paste0("../input_data_updated/",ct,"/enh_tf_list"))
  
  load(paste0("../input_data_updated/",ct,"/use_pos_list"))
  
  tf_name <- colnames(enh_tf_list[[1]])
  
  rownames(group_TF_bar) <- tf_name
  
  enriched_tf <- apply(group_TF_bar,2,function(x){names(which(x>2))})
  
  gene_group_tf.pair <- enriched_tf[ group_infor[pool$gene] ]
  
 
  pool$gene_group_index <- group_infor[pool$gene] 
  
  chr_name <- paste0("chr",c(1:22,"X"))
  
  bind_tf <- list()
  
  for(i in 1:dim(pool)[1]){
    
    tmp_enh <- enh_coord[pool$enh[i],]
    
    chr_id <- match(tmp_enh[1,1],chr_name)
    
    tmp_chr <- chr_name[chr_id]
    
    enh_pool <- subset(enh_coord,enh_coord[,1]==tmp_chr)
    
    use_id <- prodlim::row.match(tmp_enh,enh_pool)
    
    loc <- match(use_id,use_pos_list[[chr_id]])
    
    tf_count <- enh_tf_list[[chr_id]][loc,]
    
    bind_tf[[i]] <- names(which(tf_count>0))
    
  }
  
  tf_ind <- c()
  
  for(i in 1:dim(pool)[1]){
    
    tf_ind[i] <- any(bind_tf[[i]]%in%gene_group_tf.pair[[i]])
    
  }
  
  pool$tf <- mapply(function(x,y){paste0(subset(x,x%in%y),collapse = ";")},bind_tf,gene_group_tf.pair)
  
  pool$g_tf <- sapply(gene_group_tf.pair,function(x){paste0(x,collapse = ";")})
  
  
  pool$tf_ind <- tf_ind
  
  if(dim(pool)[1]>20){
  pool <- subset(pool,pool$tf_ind==1)}
  
  #check if nearest gene
  
  ng <- c()
  
  for(i in 1:dim(pool)[1]){
    
    tmp_enh <- pool[i,2]
    
    tmp_dist <- pool[i,3]
    
    all_dist <- subset(info_frame$dist,info_frame$enh==tmp_enh)
    
    if(any(all_dist<tmp_dist)){ng[i] <- 0}else{
      
      ng[i] <- 1
      
    }
    
  }
  
  pool$nearest_gene <- ng
  
  #check gene act and rank
  load("../data/gene_act_mat_updated.Rdata")
  
  pool$gene_act <- gene_act_mat[pool$gene,112]
  
  pool$gene_rank <- apply(gene_act_mat[pool$gene,],1,function(x){rank(x)[112]})
  
  #calculate degree
  
  deg <- c()
  
  pred_link <- subset(info_frame,info_frame$corrected_p>0.9)
  
  for(i in 1:dim(pool)[1]){
    
    deg[i] <- sum(pred_link$gene == pool[i,1])
    
  }
  
  pool$deg <- deg
  
  return(pool)
}





#select tier 1/3 examples: all fiters, tier 3 only supported by eqtl and Hi-C

pool <- subset(info_frame,info_frame$corrected_p>0.85&info_frame$corr>0.3&(info_frame$CaptureC==1)&info_frame$tf_count>0&info_frame$GTEX>0)
pool$gene_name <- gene_annotation[match(promoter_frame[pool[,1],1],gene_annotation[,7]),6]
pool$enh_name <- apply(enh_coord[pool$enh,],1,function(x){gsub(" ","",paste0(x[1],":",x[2],"-",x[3]))})

if(dim(pool)[1]>0){

s_pool <- select_example(pool)

tier_1 <- subset(s_pool,s_pool$gene_rank>80&s_pool$deg<11&s_pool$tf_ind==1)


tier_3 <- s_pool

}else{

tier_1 <- tier_3 <- NULL

}
#select tier 2 examples: with eqtl or Capture C

pool <- subset(info_frame,info_frame$corrected_p>0.9&info_frame$corr>0.4&(info_frame$CaptureC==1|info_frame$GTEX>0)&info_frame$tf_count>0&info_frame$dist>1.5e4)
pool$gene_name <- gene_annotation[match(promoter_frame[pool[,1],1],gene_annotation[,7]),6]
pool$enh_name <- apply(enh_coord[pool$enh,],1,function(x){gsub(" ","",paste0(x[1],":",x[2],"-",x[3]))})

tier_2 <- select_example(pool)

tier_2 <- subset(tier_2,tier_2$gene_rank>80&tier_2$deg<11)

write.table(tier_1,paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/data/examples_3_layer/auto/",ct,"_tier_1.txt"),col.names = T,row.names = F,sep="\t",quote=F)
write.table(tier_2,paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/data/examples_3_layer/auto/",ct,"_tier_2.txt"),col.names = T,row.names = F,sep="\t",quote=F)
write.table(tier_3,paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/data/examples_3_layer/auto/",ct,"_tier_3.txt"),col.names = T,row.names = F,sep="\t",quote=F)
