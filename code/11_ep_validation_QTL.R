args <- commandArgs(T)

i <- as.numeric(args[1])

eqtl_file <- args[2]

promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")
enh_frame <- read.table("../consensus_enhancer_coord.txt")
load("../data/gene_enh_frame.Rdata")
promoter_frame <- promoter_frame[,-1]

eQTL <- read.table(eqtl_file,sep="\t",header=F)
eQTL <- unique(eQTL)
#eQTL_frame <- data.frame(gene_id=eQTL$gene_id,variant_id=eQTL$variant_id)
eQTL_frame <- eQTL

generate_frame <- function(x){
  
  gene <- substr(as.character(x[1]),1,15)
  
  variant <- strsplit(as.character(x[2]),"_")[[1]]
  
  chr <- variant[1]
  
  loc <- variant[2]
  
  return(data.frame(gene=gene,chr=chr,loc=loc))
  
}

#eQTL_frame <-  apply(eQTL_frame,1,generate_frame)

#eQTL_frame <- do.call(rbind,eQTL_frame)

eQTL_frame[,1] <- as.character(eQTL_frame[,1])

eQTL_frame[,2] <- paste0("chr",as.character(eQTL_frame[,2]))

#eQTL_frame[which(eQTL_frame[,2]=="X"),2] <- 23

eQTL_frame[,3] <- as.numeric(as.character(eQTL_frame[,3]))


#load(paste0(pth,"gene_enh_frame"))
#load(paste0(pth,"enh_frame"))
#load(paste0(pth,"promoter_frame"))


gene_name <- read.table("../promoter_activity/gene_name.txt",colClasses = "character")[,1]

gene_name <- sapply(strsplit(gene_name,"[.]"),function(x){x[[1]][1]})



match <- function(x){
  
  temp_gene_name <- gene_name[as.numeric(x[1])]
  
  enh_loc <- enh_frame[as.numeric(x[2]),]
  
  ind <- which(eQTL_frame[,1]==temp_gene_name&eQTL_frame[,2]==enh_loc[,1]&eQTL_frame[,3]<enh_loc[,3]&eQTL_frame[,3]>enh_loc[,2])
  if(length(ind)==0){ind <- 0}
  return(ind)
}



n <- dim(gene_enh_frame)[1]

n_1 <- floor(n/500)


if(i < 500){
  gene_enh_frame <- gene_enh_frame[((i-1)*n_1+1):(i*n_1),]}else{
    
    gene_enh_frame <- gene_enh_frame[(499*n_1+1):n,]
    
  }

ind <- apply(gene_enh_frame,1,match)

dir_name <- gsub("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/Link_frame/eQTL/gtex","",args[2])

dir_name <- gsub(".txt","",dir_name)

opt.path <- paste0("../validation/",dir_name)

dir.create(opt.path,recursive = T)

write.table(ind,paste0(opt.path,"/",i,".txt"),col.names = F,row.names = F,quote=F)
