args <- commandArgs(TRUE)
#2 parameters:
#1: index of chunk
#2: link to the gold standards


promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")
enh_frame <- read.table("../consensus_enhancer_coord.txt")
load("../data/gene_enh_frame.Rdata")

i <- as.numeric(args[1])
# 
# enh_center <- (enh_frame[,2]+enh_frame[,3])/2
# 
# enh_frame[,2] = enh_center-1000
# 
# enh_frame[,3] = enh_center+1000

#promoter_frame[,1] <- paste0("chr",promoter_frame[,1])

promoter_frame <- promoter_frame[,-1]

n <- dim(gene_enh_frame)[1]

n_1 <- floor(n/500)

match <- function(x){
  
  gene_loc <- promoter_frame[as.numeric(x[1]),]
  
  enh_loc <- enh_frame[as.numeric(x[2]),]
  
  ind<- which(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&link_frame[,1]==gene_loc[,1]&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))
  
  if(length(ind)>0){
  return(paste0(ind,collapse = ";"))}else{
    
    return(0)
    
  }
  
}



link_frame <- read.table(args[2])

link_frame[,1] <-as.character(link_frame[,1])

link_frame[,1] <- paste0("chr",link_frame[,1])

#library(parallel)

if(i < 500){
  gene_enh_frame <- gene_enh_frame[((i-1)*n_1+1):(i*n_1),]}else{
    
    gene_enh_frame <- gene_enh_frame[(499*n_1+1):n,]
    
  }

ind <- apply(gene_enh_frame,1,match)

dir_name <- gsub("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/data/Link_frame/hQTL/","",args[2])
dir_name <- gsub(".txt","",dir_name)

opt.path <- paste0("../validation/",dir_name)

dir.create(opt.path,recursive = T)

write.table(ind,paste0(opt.path,"/",i,".txt"),col.names = F,row.names = F,quote=F)
