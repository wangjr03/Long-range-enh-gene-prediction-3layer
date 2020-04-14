args <- commandArgs(TRUE)


pth <- paste0("../IM-PET/",args[3],"/result/all_pairs.txt")

gene_enh_frame <- data.table::fread(pth)
gene_enh_frame <- as.data.frame(gene_enh_frame)
i <- as.numeric(args[1])

n <- dim(gene_enh_frame)[1]

n_1 <- floor(n/500)

match <- function(x){
  
  gene_loc <- x[c(6,7)]
  
  enh_loc <- x[1:3]
  
  ind<- any(link_frame[,1]==as.character(enh_loc[1])&link_frame[,2]<as.numeric(enh_loc[3])&link_frame[,3]>as.numeric(enh_loc[2])&
              link_frame[,1]==as.character(gene_loc[1])&link_frame[,4]<as.numeric(gene_loc[2])&link_frame[,5]>as.numeric(gene_loc[2]))+any(link_frame[,1]==as.character(enh_loc[1])&link_frame[,4]<as.numeric(enh_loc[3])&link_frame[,5]>as.numeric(enh_loc[2])&
                                                                                                                                             link_frame[,1]==gene_loc[1]&link_frame[,2]<as.numeric(gene_loc[2])&link_frame[,3]>as.numeric(gene_loc[2]))
  
  
  return(ind)
  
}



link_frame <- read.table(args[2])

link_frame[,1] <-as.character(link_frame[,1])

library(parallel)

if(i < 500){
  gene_enh_frame <- gene_enh_frame[((i-1)*n_1+1):(i*n_1),]}else{
    
    gene_enh_frame <- gene_enh_frame[(499*n_1+1):n,]
    
  }

gene_enh_frame <- as.data.frame(gene_enh_frame)

ind <- apply(gene_enh_frame,1,match)

dir_name <- gsub("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/Link_frame/","",args[2])
dir_name <- gsub(".txt","",dir_name)

o.pth <- paste0("../IM-PET/",args[3],"/",dir_name)

dir.create(o.pth,recursive=T)

write.table(ind,paste0(o.pth,"/",i,".txt"),col.names = F,row.names = F,quote=F)
