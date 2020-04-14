exons <- read.table("../gene_annotation/exons_gencode_v19_protein_coding_gene.txt",header=T)

gene_name <- read.table("../promoter_activity/gene_name.txt")[,1]

id <- sapply(gene_name,function(x){
  
  which(as.character(exons[,4])==x)
  
})

gene_exp_mat <- list()



for(i in 1:127){
  print(i)
  file <- paste0("../gene_expr_exons/",i,"_RNA_seq.Rdata")
  
  load(file)
  
  #exons$score <- average_score
  
  score <- list()

  for(j in 1:23){
  
  score <- c(score,score_chr[[j]])
  
  } 

  gene_exp <- sapply(id,function(x){
    
    #if(length(x)==0){return(0)}else{
    
    
    
      tmp_score <- do.call(rbind,score[x])
      if(dim(tmp_score)[2]>0){      
      return(mean(tmp_score[,2]))}else{return(0)}
    
  })
  
  gene_exp_mat[[i]] <- gene_exp
  
}


gene_exp_mat <- do.call(rbind,gene_exp_mat)
gene_exp_mat <- as.data.frame(t(gene_exp_mat))

save(gene_exp_mat,file="../data/gene_act_mat_updated.Rdata")



