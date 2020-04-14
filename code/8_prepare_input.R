#pre process step
#generate enh tf matrix
#2 parameters
#1: cell type
#2: max distance between enhancer and promoter
#3: correlation version, promoter or exons

args <- commandArgs(T)

cell_type <- as.numeric(args[1])

output_path <- paste0("../input_data_",args[3],"/",cell_type,"/")

distance <- as.numeric(args[2])

dir.create(output_path,recursive = T)

setwd(output_path)

#so far should have a matrix for each chr containing TF counts of each enhancer
#now get use pos of each chr
#filter out inactive TFs, also update use-pos
load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/all_tf_name")
tf_name <- unique(all_tf_name)

load(paste0("../../data/gene_act_mat_",args[3],".Rdata"))

tf_exp <- gene_act_mat

ensg_name <- read.table("../../promoter_activity/gene_name.txt")

gene_annotation <- read.table("../../gene_annotation/gene_annotation_V19.txt")

tf_ensg_id <- as.character(gene_annotation[match(as.character(tf_name),as.character(gene_annotation[,6])),7])

tf_loc <- match(tf_ensg_id,as.character(ensg_name[,1]))

expressed_tf <- subset(tf_name,floor(tf_exp[tf_loc,cell_type])>0)

tf_name <- expressed_tf



print("generating enhancer-TF matrix...")
use_pos_list <- list()

enh_tf_list <- list()

for(i in 1:23){
  
  if(i!=23){
    
    load(paste0("../../enh_tf_Rdata/chr",i,".Rdata"))
    
    enh_tf <- enh_tf[,which(colnames(enh_tf)%in%tf_name)]
    
    zero_index <- which(apply(enh_tf,1,sum)==0)
    
    if(length(zero_index)>0){
      
      enh_tf <- enh_tf[-zero_index,]
      
      use_pos <- use_pos[-zero_index]
      
    }
   
    use_pos_list[[i]] <- use_pos
    enh_tf_list[[i]] <- enh_tf
    
  }else{
    
    load(paste0("../../enh_tf_Rdata/chr","X",".Rdata"))
    
    enh_tf <- enh_tf[,which(colnames(enh_tf)%in%tf_name)]
    
    zero_index <- which(apply(enh_tf,1,sum)==0)
    
    if(length(zero_index)>0){
      
      enh_tf <- enh_tf[-zero_index,]
      
      use_pos <- use_pos[-zero_index]
      
    }
    
    use_pos_list[[i]] <- use_pos
    enh_tf_list[[i]] <- enh_tf
    
  }

}



save(enh_tf_list,file="enh_tf_list")

save(use_pos_list,file="use_pos_list")
print("done")
#load enh_tf list
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/exp_enh_tf_list/CD4+.Rdata")

#load pot_dist_list

#print("generate potential e-p pairs...")
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_act_list")
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_corr_list_qt_spearman")
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_enh_list")
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_dist_list")
# 
# 
# tf_pot_dist_list <- tf_pot_corr_list <- tf_pot_enh_list <-tf_pot_act_list <- list()
# 
# for(i in 1:23){
#   
#   empt_gene <- c()
#   
#   tf_pot_dist_list[[i]]<-tf_pot_act_list[[i]]<- tf_pot_corr_list[[i]] <- tf_pot_enh_list[[i]] <- list()
# 
#   for(j in 1:length(pot_enh_list[[i]])){
# 
#     index <- which(pot_enh_list[[i]][[j]]%in%use_pos_list[[i]])
#     
#     if(length(index)==0){empt_gene <- c(empt_gene,j)}
#     
#     tf_pot_enh_list[[i]][[j]] <- pot_enh_list[[i]][[j]][index]
# 
#     tf_pot_corr_list[[i]][[j]] <- pot_corr_list[[i]][[j]][index]
# 
#     tf_pot_dist_list[[i]][[j]] <- pot_dist_list[[i]][[j]][index]
# 
#     tf_pot_act_list[[i]][[j]] <- pot_act_list[[i]][[j]][index]
#   }
# 
#   names(tf_pot_enh_list[[i]]) <- names(tf_pot_dist_list[[i]]) <- names(tf_pot_corr_list[[i]]) <- names(pot_dist_list[[i]])
# 
#   if(length(empt_gene)>0){
#   
#   tf_pot_act_list[[i]] <- tf_pot_act_list[[i]][-empt_gene]
#   
#   tf_pot_corr_list[[i]] <- tf_pot_corr_list[[i]][-empt_gene]
#   
#   tf_pot_dist_list[[i]] <- tf_pot_dist_list[[i]][-empt_gene]
#   
#   tf_pot_enh_list[[i]] <- tf_pot_enh_list[[i]][-empt_gene]
#   }
# }
# 
# 
# pot_act_list <- tf_pot_act_list
# 
# pot_enh_list <- tf_pot_enh_list
# 
# pot_dist_list <- tf_pot_dist_list
# 
# pot_corr_list <- tf_pot_corr_list
# 
# save(pot_enh_list,file=paste0(output_path,"pot_enh_list"))
# 
# save(pot_act_list,file=paste0(output_path,"pot_act_list"))
# 
# save(pot_dist_list,file=paste0(output_path,"pot_dist_list"))
# 
# save(pot_corr_list,file=paste0(output_path,"pot_corr_list"))

print(paste("generating potential e-p pairs using window size +/-",distance))

#load("../../data/pot_act_list")
load("../../data/pot_dist_list.Rdata")
load("../../data/pot_enh_list.Rdata")
load(paste0("../../data/pot_corr_list_",args[3],".Rdata"))

get_ind <- function(x,y){
  
  ind <- which(!y>distance)
  
  return(x[ind])
}

for(i in 1:23){
  
  pot_enh_list[[i]] <- mapply(get_ind,pot_enh_list[[i]],pot_dist_list[[i]])
  #pot_act_list[[i]] <- mapply(get_ind,pot_act_list[[i]],pot_dist_list[[i]])
  pot_corr_list[[i]] <- mapply(get_ind,pot_corr_list[[i]],pot_dist_list[[i]])
  pot_dist_list[[i]] <- mapply(get_ind,pot_dist_list[[i]],pot_dist_list[[i]])
}

#deal with empty gene


for(i in 1:23){
  
  id <- which(sapply(pot_enh_list[[i]],length)!=0)
  
  pot_enh_list[[i]] <- pot_enh_list[[i]][id]
  
  #pot_act_list[[i]] <- pot_act_list[[i]][id]
  
  pot_corr_list[[i]] <- pot_corr_list[[i]][id]
  
  pot_dist_list[[i]] <- pot_dist_list[[i]][id]
  
  
}




#save(pot_act_list,file="pot_act_list")
save(pot_enh_list,file="pot_enh_list")
save(pot_corr_list,file="pot_corr_list")
save(pot_dist_list,file="pot_dist_list")


#system(paste("cp /mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_act_list",paste0(output_path,"pot_act_list")))
#system(paste("cp /mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_enh_list",paste0(output_path,"pot_enh_list")))
#system(paste("cp /mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_dist_list",paste0(output_path,"pot_dist_list")))
#system(paste("cp /mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/pot_corr_list_qt_spearman",paste0(output_path,"pot_corr_list")))


print("done")

print("generate tf count data for each enhancer")
#generate tf count list
load("pot_enh_list")
tf_count_list <- list()

enh_tf_count <- sapply(enh_tf_list,function(x){apply(x,1,sum)})

for(i in 1:23){
  
  tf_count_list[[i]] <- list()
  
  for(j in 1:length(pot_enh_list[[i]])){
    
    tf_count_list[[i]][[j]] <- enh_tf_count[[i]][match(pot_enh_list[[i]][[j]],use_pos_list[[i]])]
    tf_count_list[[i]][[j]][which(is.na(tf_count_list[[i]][[j]]))] <- 0
    
  }
  
}



save(tf_count_list,file="tf_count_list")



#gene_act
print("generate activity data...")

gene_infor_matrix <- c()


for(i in 1:23){
  
  gene_id <- names(pot_enh_list[[i]])
  
  count <- length(pot_enh_list[[i]])
  
  tmp <- data.frame(gene_id=gene_id,chr=i,id=1:count)
  
  gene_infor_matrix <- rbind(gene_infor_matrix,tmp)
}


gene_name_list <- sapply(pot_enh_list,function(x){names(x)})

save(gene_infor_matrix,file="gene_infor_matrix")

save(gene_name_list,file="gene_name_list")


all_gene_name <- read.table("../../promoter_activity/gene_name.txt")

all_gene_exp <- gene_act_mat

gene_loc <- match(gene_infor_matrix[,1],all_gene_name[,1])

gene_act <- floor(all_gene_exp[gene_loc,cell_type])

gene_act <- sapply(split(gene_act,gene_infor_matrix[,2]),as.list)


#enh_act

all_enhancer_name <- read.table("../../consensus_enhancer_coord.txt")

load("../../data/enh_act_mat.Rdata")

all_enhancer_act <- enh_act_mat

enh_act_list <- list()

for(i in 1:23){
  
  enh_act_list[[i]] <- list()
  
  if(i !=23){
    
    temp_enh <- subset(all_enhancer_name,all_enhancer_name[,1]==paste0("chr",i))
    
  }else{
    
    temp_enh <- subset(all_enhancer_name,all_enhancer_name[,1]=="chrX")
  }
  
  temp_loc <- prodlim::row.match(temp_enh,all_enhancer_name)
  
  temp_act <- floor(all_enhancer_act[temp_loc,cell_type])
  
  enh_act_list[[i]] <- sapply(pot_enh_list[[i]],function(x){as.vector(temp_act[x])})
  
  names(enh_act_list[[i]]) <- names(pot_enh_list[[i]])
  
}


save(gene_act,file = "gene_act")

save(enh_act_list,file="enh_act_list")

pot_act_list <- enh_act_list

save(pot_act_list,file="pot_act_list")

print("done")


#generate enh_act_chr
enh_act_chr <- list()


for(i in 1:23){
  
  enh_act_chr[[i]] <- c()
  
  if(i !=23){
    
    temp_enh <- subset(all_enhancer_name,all_enhancer_name[,1]==paste0("chr",i))
    
  }else{
    
    temp_enh <- subset(all_enhancer_name,all_enhancer_name[,1]=="chrX")
  }
  temp_loc <- prodlim::row.match(temp_enh,all_enhancer_name)
  
  temp_act <- floor(all_enhancer_act[temp_loc,cell_type])
  
  enh_act_chr[[i]] <- temp_act
  
}

save(enh_act_chr,file="enh_act_chr")

#generate gene_act_chr

gene_act_chr <- sapply(gene_act,unlist)

save(gene_act_chr,file="gene_act_chr")


# #for each gene, get active pair
# act_pot_dist_list <- act_pot_corr_list <- act_pot_enh_list <- list()
# 
# for(i in 1:23){
#   
#   act_pot_dist_list[[i]] <- act_pot_corr_list[[i]] <- act_pot_enh_list[[i]] <- list()
#   
#   for(j in 1:length(pot_enh_list[[i]])){
#     
#     index <- which(pot_act_list[[i]][[j]]==1)
#     
#     act_pot_enh_list[[i]][[j]] <- pot_enh_list[[i]][[j]][index]
#     
#     act_pot_corr_list[[i]][[j]] <- pot_corr_list[[i]][[j]][index]
#     
#     act_pot_dist_list[[i]][[j]] <- pot_dist_list[[i]][[j]][index]
#     
#   }
#   
#   names(act_pot_enh_list[[i]]) <- names(act_pot_dist_list[[i]]) <- names(act_pot_corr_list[[i]]) <- names(pot_dist_list[[i]])
#   
# }
# 
# save(act_pot_enh_list,file="/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/activity_filter_data/act_pot_enh_list")
# 
# save(act_pot_dist_list,file="/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/activity_filter_data/act_pot_dist_list")
# 
# save(act_pot_corr_list,file="/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/activity_filter_data/act_pot_corr_list")


#initialize interaction based on correlations
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/corr_dist_jr_predict")
# 
# inter_d <- density(unlist(group_corr_data))
# 
# non_inter_d <- density(unlist(group_non_corr_data))
# 
# inter_prob <- function(tmp.p){
#   
#   loc <- which(inter_d$x>tmp.p)[1]
#   
#   if(is.na(loc)){
#     
#     return(0)
#     
#   }else if(loc>1){
#     
#     return(inter_d$y[loc-1])
#     
#   }else if(loc==1){
#     
#     return(0)
#     
#   }
#   
# }
# 
# non_inter_prob <- function(tmp.p){
#   
#   loc <- which(non_inter_d$x>tmp.p)[1]
#   
#   if(is.na(loc)){
#     
#     return(0)
#     
#   }else if(loc>1){
#     
#     return(non_inter_d$y[loc-1])
#     
#   }else if(loc==1){
#     
#     return(0)
#     
#   }
#   
# }
# 
# 
# corr_init <- function(x){
#   
#   p_i <- inter_prob(x)
#   
#   p_n <- non_inter_prob(x)
#   
#   p <- p_i/(p_i+p_n+1e-10)
#   
#   index <- c(0,1)p[which.max(1-p,p)]
#   
#   return(index)
# }
# 
# 
# corr_init <- Vectorize(corr_init)
# 
# inter_prob <- Vectorize(inter_prob)
# 
# non_inter_prob <- Vectorize(non_inter_prob)

#sample index for each active pair

print("initialize e-p link...")

load("pot_corr_list")
load("pot_dist_list")
load("pot_enh_list")

corr_init <- function(x){
  
  x <- max(x,0)
  
  p_l <- x
  p_n <- 1-x
  
  index <- sample(c(0,1),1,prob = c(p_n,p_l))
  return(index)
}

corr_init <- Vectorize(corr_init)

index_list <- lapply(pot_corr_list,function(x){
  
  lapply(x,corr_init)
  
})

inter_enh_list <- list()

inter_corr_list <- list()

inter_dist_list <- c()

for(i in 1:23){
  
  inter_enh_list[[i]] <- list()
  
  inter_dist_list[[i]] <- list()
  
  inter_corr_list[[i]] <- list()
  
  for(j in 1:length(pot_enh_list[[i]])){
    
    inter_enh_list[[i]][[j]] <- pot_enh_list[[i]][[j]][which(index_list[[i]][[j]]==1&tf_count_list[[i]][[j]]>0&(gene_act[[i]][[j]]*pot_act_list[[i]][[j]])>0)]
    
    inter_corr_list[[i]][[j]] <- pot_corr_list[[i]][[j]][which(index_list[[i]][[j]]==1&tf_count_list[[i]][[j]]>0&(gene_act[[i]][[j]]*pot_act_list[[i]][[j]])>0)]
    
    inter_dist_list[[i]][[j]] <- pot_dist_list[[i]][[j]][which(index_list[[i]][[j]]==1&tf_count_list[[i]][[j]]>0&(gene_act[[i]][[j]]*pot_act_list[[i]][[j]])>0)]
    
  }
}

save(inter_enh_list,file="inter_enh_list")
save(inter_corr_list,file="inter_corr_list")

print("done")
#generate TF_GE matrix
print("generate gene-TF matrix...")
load("enh_tf_list")

load("use_pos_list")

load("inter_enh_list")

load("pot_enh_list")

TF_GE_list <- list()

for(i in 1:23){
  
  tmp_TF_GE <- c()
  
  for(j in 1:length(pot_enh_list[[i]])){
    
    inter.enh <- inter_enh_list[[i]][[j]]
    
    inter.enh.loc <- match(inter.enh,use_pos_list[[i]])
    
    
    if(length(inter.enh.loc)==0){
      inter.tf.mat <- enh_tf_list[[i]][1,]
    }else{
      inter.tf.mat <- enh_tf_list[[i]][inter.enh.loc,]}
    
    if(class(inter.tf.mat)!="numeric"){
      gene.tf <- apply(inter.tf.mat,2,mean)}else{
        gene.tf <- inter.tf.mat
      }
    
    tmp_TF_GE <- rbind(tmp_TF_GE,gene.tf)
    #print(j)
  }
  
  rownames(tmp_TF_GE) <- names(pot_enh_list[[i]])
  
  TF_GE_list[[i]] <- tmp_TF_GE
  
}

all_TF_GE <- do.call(rbind,TF_GE_list)

gene_ind <- which(apply(all_TF_GE,1,sum)==0)

mean_tf <- apply(all_TF_GE,2,mean)

for(i in gene_ind){
  
  all_TF_GE[i,] <- mean_tf
  
}

TF_GE_list_ind <- lapply(TF_GE_list,function(x){
  ind <- which(apply(x,1,sum)==0)
  return(ind)
})

for(i in 1:23){
  
  if(length(TF_GE_list_ind[[i]])==0){next}
  for(j in TF_GE_list_ind[[i]]){
    
    TF_GE_list[[i]][j,] <- mean_tf
    
  }
  
  
}

save(all_TF_GE,file="all_TF_GE")

save(TF_GE_list,file="TF_GE_list")
print("done")
#having TF_GE_matrix, generate a hierarcical clustering tree
library(pheatmap)
ind <- which(apply(all_TF_GE,2,sum)>0)

gene_ind <- which(apply(all_TF_GE,1,sum)==0)

all_TF_GE[gene_ind,] <- apply(all_TF_GE,2,mean)

print("generate gene group")
phm <- pheatmap(all_TF_GE[,ind],color = colorRampPalette(c("snow", "red"))(50),clustering_distance_row = "correlation", 
                clustering_distance_col = "correlation", fontsize=9, fontsize_row=1, fontsize_col=3,silent=T)


tree_gene <- phm$tree_row

save(tree_gene,file="tree_gene")


system("cp /mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/No_conservation_motif/all_inter_dist ./")


print("finish!")







