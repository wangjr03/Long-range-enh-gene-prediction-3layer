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
  
  exons$score <- average_score
  
  gene_exp <- sapply(id,function(x){
    
    if(length(x)==0){return(0)}else{
    
    mean(exons[x,6])}
    
  })
  
  gene_exp_mat[[i]] <- gene_exp
  
}


gene_exp_mat <- as.data.frame(gene_exp_mat)

save(gene_exp_mat,file="../data/gene_act_mat_exons.Rdata")



