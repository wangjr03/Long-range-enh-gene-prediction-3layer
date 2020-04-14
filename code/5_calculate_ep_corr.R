library(preprocessCore)
#1 parameters
#1:specify promoter or exons, correspond to rna-seq version

args <- commandArgs(T)


calculate_cor <- function(x,y){
  
  p <- cor(as.numeric(tmp_gene_act[x,]),as.numeric(tmp_enh_act[y,]))
  
  return(p)
  
}

calculate_cor <- Vectorize(calculate_cor,vectorize.args = "y")


#system("cp ../promoter_activity/RNA_seq_127_matrix.txt ../data/")
load(paste0("../data/gene_act_mat_",args[1],".Rdata"))


#load enhancer/promoter activity profile
load("../data/enh_act_mat.Rdata")
enh_coord <- read.table("../consensus_enhancer_coord.txt")


#gene_act_mat <- read.table("../data/RNA_seq_127_matrix.txt")

promoter <- read.table("../promoter_activity/promoter_frame.txt")

#promoter <- read.table("../promoter_activity/promoter_frame.txt")

promoter[,2] <- as.character(promoter[,2])

promoter[which(promoter[,2]=="chr23"),2] <- "chrX"

names_level <- paste0("chr",c(1:22,"X"))

enh_coord[,1] <- factor(enh_coord[,1],levels = names_level)

promoter[,2] <- factor(promoter[,2],levels = names_level)

#quantile normalize activities
enh_act_mat <- normalize.quantiles(as.matrix(enh_act_mat))
gene_act_mat <- normalize.quantiles(as.matrix(gene_act_mat))

a <- apply(gene_act_mat,1,sd)

a  <- which(a==0)

gene_act_mat[a,] <- gene_act_mat[a,]+abs(rnorm(length(a)*127,0,1e-5))

#get act profile in each chromosome
enh_act_chr <- gene_act_chr <- list()

k <- 1

for(i in names_level){
  
  enh_id <- which(enh_coord[,1]==i)
  
  gene_id <- which(promoter[,2]==i)
  
  enh_act_chr[[k]] <- enh_act_mat[enh_id,]
  
  gene_act_chr[[k]] <- gene_act_mat[gene_id,]
  
  k <- k+1
  
}


#load all pot enhancer promoter pool
load("../data/pot_enh_list.Rdata")

pot_corr_list <- list()

for(i in 1:23){
  
  tmp_enh_act <- enh_act_chr[[i]]
  
  tmp_gene_act <- gene_act_chr[[i]]
  
  tmp_pot_enh <- pot_enh_list[[i]]
  
  n <- length(tmp_pot_enh)
  
  pot_corr_list[[i]] <- sapply(1:n,function(x){
    
    calculate_cor(x,tmp_pot_enh[[x]])
    
  })
  
  
}

save(pot_corr_list,file=paste0("../data/pot_corr_list_",args[1],".Rdata"))











