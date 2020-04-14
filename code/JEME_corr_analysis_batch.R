args <- commandArgs(T)

ct <- as.numeric(args[1])


library(ppcor)

load(paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/final_prediction/",ct,"/info_frame.Rdata"))

pred <- subset(info_frame,info_frame$corrected_p>0.9)



#calculate correlation and partial correlation for JEME

jeme_pair <- read.csv(paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/Yip/prediction/encoderoadmap_lasso/encoderoadmap_lasso.",ct,".csv"),header=F)

jeme_enh <- as.data.frame(do.call(rbind,strsplit(as.character(jeme_pair[,1]),"[:]|-")))

dir.create(paste0("../mapping/",ct),recursive=T)

setwd(paste0("../mapping/",ct))

write.table(jeme_enh,"jeme_enh.txt",col.names = F,row.names = F,sep="\t",quote=F)

system("python /mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/Yip/intermediate_result/0929/overlap_in_two_files2.py -i1 jeme_enh.txt -i2 /mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/consensus_enhancer_coord.txt -o1 jeme_enh_mapping.txt")


file.1 <- "./jeme_enh_mapping.txt"

con.1 <- file(file.1,"r")

all_line <- readLines(con.1)

whole_mapping <- list()

jeme_unique_enh <- c()

for(i in 1:length(all_line)){
  
  print(i)
  line <- strsplit(all_line[i],"\t")
  if ( length(line) == 0 ) {
    break
  }
  
  temp_jeme_enh <- line[[1]][1:3]
  
  jeme_unique_enh <- rbind(jeme_unique_enh,temp_jeme_enh)
  
  n <- length(line[[1]])
  
  for(j in 4:n){
    
    my_enh <- line[[1]][j]
    
    my_enh <- gsub(" ","",my_enh)
    
    my_enh <- strsplit(my_enh,"[(]|[)]|,|[']")[[1]][c(3,5,6)]
    
    temp_mapping <- c(temp_jeme_enh,my_enh,i)
    
    whole_mapping[[i]] <- temp_mapping
  }
  
  
}

whole_mapping <- do.call(rbind,whole_mapping)

whole_mapping <- as.data.frame(whole_mapping)

jeme_pair <- cbind(jeme_enh,jeme_pair[,c(2,3)])

in_id <- which(!is.na(prodlim::row.match(jeme_pair[,c(1,2,3)],whole_mapping[,c(1,2,3)])))

jeme_intersect <- jeme_pair[in_id,]


enh_frame <- read.table("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj9_SNP_network/data/consensus_enhancer_coord.txt")
promoter_frame <- read.table("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj9_SNP_network/data/promoter.txt")
gene_anno <- read.table("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj9_SNP_network/data/gene_name.txt")

#get intersect gene name
jeme_gene_name <- unique(sapply(strsplit(as.character(jeme_intersect[,4]),"[$]"),function(x){x[1]}))

#jeme_intersect[,4] <- sapply(strsplit(as.character(jeme_intersect[,4]),"[.]"),function(x){x[1]})
jeme_intersect[,4] <- sapply(strsplit(as.character(jeme_intersect[,4]),"[$]"),function(x){x[1]})


common_gene <- intersect(jeme_gene_name,gene_anno[,1])

#get jeme common pair

jeme_common <- subset(jeme_intersect,jeme_intersect[,4]%in%common_gene)


#expand jeme pair: for some enhancer do replecation

whole_mapping <- as.data.frame(whole_mapping)

whole_mapping[,1] <- as.character(whole_mapping[,1])

whole_mapping[,2] <- as.numeric(as.character(whole_mapping[,2]))

whole_mapping[,3] <- as.numeric(as.character(whole_mapping[,3]))

whole_mapping[,4] <- as.character(whole_mapping[,4])

whole_mapping[,5] <- as.numeric(as.character(whole_mapping[,5]))

whole_mapping[,6] <- as.numeric(as.character(whole_mapping[,6]))


jeme_common[,1] <- as.character(jeme_common[,1])

jeme_common[,2] <- as.numeric(as.character(jeme_common[,2]))

jeme_common[,3] <- as.numeric(as.character(jeme_common[,3]))

jeme_expnd <- list()

for(i in 1:dim(jeme_common)[1]){
  if(i%%100==0){print(i)}
  tmp_enh <- jeme_common[i,c(1,2,3)]
  
  mapped_id <- which(whole_mapping[,1]==tmp_enh[,1]&whole_mapping[,2]==tmp_enh[,2]&whole_mapping[,3]==tmp_enh[,3])
  
  tmp_frame <- data.frame(gene=jeme_common[i,4],enh=whole_mapping[mapped_id,c(4,5,6)],score = jeme_common[i,5])
  
  jeme_expnd[[i]] <-tmp_frame
  
}

jeme_expnd <- do.call(rbind,jeme_expnd)

jeme_expnd$enh_id <- prodlim::row.match(jeme_expnd[,c(2,3,4)],enh_frame)

jeme_expnd$gene_id <- match(jeme_expnd[,1],as.character(promoter_frame[,4]))

jeme_expnd <- subset(jeme_expnd,jeme_expnd$score>0.3)

jeme_expnd$gene <- as.character(jeme_expnd$gene)

jeme_pool <- split(jeme_expnd,jeme_expnd$gene)

jeme_enh_pair <- sapply(jeme_pool,function(x){
  
  if(length(x$enh_id)>1){
    
    tmp <- as.data.frame(t(combn(x$enh_id,2)))
    tmp$gene <- unique(x$gene_id)
    return(tmp)
  }
  
  
  
})

jeme_enh_pair <- do.call(rbind,jeme_enh_pair)

load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj9_SNP_network/data/gene_act_mat_updated.Rdata")
load("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj9_SNP_network/data/enh_act_mat.Rdata")

jeme_cor <- apply(jeme_enh_pair,1,function(x){
  
  cor(as.numeric(enh_act_mat[x[1],]),as.numeric(enh_act_mat[x[2],]))
  
})








jeme_pcor <- apply(jeme_enh_pair,1,function(x){
  
  p_frame <- data.frame(as.numeric(enh_act_mat[as.numeric(x[1]),]),as.numeric(enh_act_mat[as.numeric(x[2]),]),as.numeric(gene_act_mat[as.numeric(x[3]),]))
  
  tryCatch(pcor(p_frame,method="spearman")$estimate[1,2],error=function(e){return(logical(0))})
  
})

jeme_pcor <- unlist(jeme_pcor)




my_pool <- split(pred,pred$gene)

my_enh_pair <- sapply(my_pool,function(x){
  
  if(length(x$enh)>1){
    
    tmp <- as.data.frame(t(combn(x$enh,2)))
    tmp$gene <- unique(x$gene)
    return(tmp)
  }
  
  
  
})

my_enh_pair <- do.call(rbind,my_enh_pair)

my_cor <- apply(my_enh_pair,1,function(x){
  
  cor(as.numeric(enh_act_mat[x[1],]),as.numeric(enh_act_mat[x[2],]))
  
})



my_pcor <- apply(my_enh_pair[sample(1:dim(my_enh_pair)[1],length(jeme_pcor)),],1,function(x){
  
  p_frame <- data.frame(as.numeric(enh_act_mat[as.numeric(x[1]),]),as.numeric(enh_act_mat[as.numeric(x[2]),]),as.numeric(gene_act_mat[as.numeric(x[3]),]))
  
  tryCatch(pcor(p_frame,method="spearman")$estimate[1,2],error=function(e){return(logical(0))})
  
})

my_pcor <- unlist(my_pcor)



rd_pair <- my_enh_pair

rd_pair[,1] <- sample(rd_pair[,1])

rd_pair[,2] <- sample(rd_pair[,2])

rd_pair[,3] <- sample(rd_pair[,3])

rd_cor <- apply(rd_pair,1,function(x){
  
  cor(as.numeric(enh_act_mat[x[1],]),as.numeric(enh_act_mat[x[2],]))
  
})

rd_pcor <- apply(rd_pair[sample(1:dim(my_enh_pair)[1],length(jeme_pcor)),],1,function(x){
  
  p_frame <- data.frame(as.numeric(enh_act_mat[as.numeric(x[1]),]),as.numeric(enh_act_mat[as.numeric(x[2]),]),as.numeric(gene_act_mat[as.numeric(x[3]),]))
  
  pcor(p_frame)$estimate[1,2]
  
})

dir.create(paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/output/correlation_analysis/",ct))

setwd(paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/output/correlation_analysis/",ct))

save(rd_cor,jeme_cor,my_cor,file="correlation.Rdata")

save(rd_pcor,jeme_pcor,my_pcor,file="partial_correlation.Rdata")



