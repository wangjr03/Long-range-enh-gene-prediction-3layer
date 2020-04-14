setwd("../gene_expr_promoter/")

file <- dir()

gene_act_mat <- list()

for(i in 1:127){
  
  load(paste0(i,"_RNA_seq.Rdata"))
  
  gene_act_mat[[i]] <- average_score
  
}

gene_act_mat <- do.call(rbind,gene_act_mat)

gene_act_mat <- as.data.frame(gene_act_mat)

gene_act_mat <- t(gene_act_mat)

save(gene_act_mat,file="../data/gene_act_mat_promoter.Rdata")
