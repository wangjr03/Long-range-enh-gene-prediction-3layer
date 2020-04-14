args <- commandArgs(T)

ct <- as.numeric(args[1])


#define gold standard
load(paste0("../output_exons_updated/",ct,"/info_frame.Rdata"))

gs_pth <- args[2]

gs_name <- args[3]

setwd(gs_pth)

n <- length(dir())

ind <- list()
for(i in 1:n){
  
  ind[[i]] <- data.table::fread(paste0(i,".txt"))
  
}
ind <- do.call(rbind,ind)
ind <- as.data.frame(ind)

info_frame$CaptureC <- ind[,1]


setwd("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/code")


jeme_pair <- read.csv(paste0("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/Yip/JEME-master/JEME_code/Roadmap/2_second_step_modeling/result/",ct,".all.sw.se.sp.csv.final.arff.apply.predictions.csv.out"),header=F)

jeme_enh <- as.data.frame(do.call(rbind,strsplit(as.character(jeme_pair[,1]),"[:]|-")))

dir.create(paste0("../../compare_JEME/",ct),recursive=T)

setwd(paste0("../../compare_JEME/",ct))

write.table(jeme_enh,paste0("../../compare_JEME/",ct,"/jeme_enh.txt"),col.names = F,row.names = F,sep="\t",quote=F)

cmd <- "python /mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/Yip/intermediate_result/0929/overlap_in_two_files2.py -i1 jeme_enh.txt -i2 ../../ENCODE_ROADMAP_127/consensus_enhancer_coord.txt -o1 jeme_enh_mapping.txt"

system(cmd)

setwd("../../ENCODE_ROADMAP_127/code/")

file.1 <- paste0("../../compare_JEME/",ct,"/jeme_enh_mapping.txt")

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


#get JEME pairs that are related with these pairs
jeme_pair <- cbind(jeme_enh,jeme_pair[,c(2,3)])



in_id <- which(!is.na(prodlim::row.match(jeme_pair[,c(1,2,3)],whole_mapping[,c(1,2,3)])))


jeme_intersect <- jeme_pair[in_id,]


#get my enhancer set

enh_frame <- read.table("../consensus_enhancer_coord.txt")
promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")
gene_anno <- read.table("../promoter_activity/gene_name.txt")

#get intersect gene name
jeme_gene_name <- unique(sapply(strsplit(as.character(jeme_intersect[,4]),"[$]"),function(x){x[1]}))

#jeme_intersect[,4] <- sapply(strsplit(as.character(jeme_intersect[,4]),"[.]"),function(x){x[1]})
jeme_intersect[,4] <- sapply(strsplit(as.character(jeme_intersect[,4]),"[$]"),function(x){x[1]})


common_gene <- intersect(jeme_gene_name,gene_anno[,1])

#get jeme common pair

jeme_common <- subset(jeme_intersect,jeme_intersect[,4]%in%common_gene)


#get my common pair

common_gene_id <- match(common_gene,gene_anno[,1])

common_enh_id <- which(!is.na(prodlim::row.match(enh_frame,unique(whole_mapping[,c(4,5,6)]))))

my_common_pair <- subset(info_frame,info_frame$gene%in%common_gene_id&info_frame$enh%in%common_enh_id&info_frame$dist<1.1e6)


tmp_data <- list()

for(i in 1:dim(my_common_pair)[1]){
  if(i%%100==0){print(i)}
  gene_id <- as.character(gene_anno[my_common_pair[i,1],1])
  
  enh_loc <- enh_frame[my_common_pair[i,2],]
  
  tmp_data[[i]] <- data.frame(gene_id,enh_loc)
  
  
}

tmp_data <- do.call(rbind,tmp_data)

my_common_pair <- cbind(my_common_pair,tmp_data)

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

#find common pair

names(jeme_expnd) <- c("gene_name","enh_chr","enh_start","enh_stop","jeme_score")

names(my_common_pair)[(dim(my_common_pair)[2]-3):(dim(my_common_pair)[2])] <- c("gene_name","enh_chr","enh_start","enh_stop")

#remove proximal pair
common_pair <- merge(jeme_expnd,my_common_pair,by=c("gene_name","enh_chr","enh_start","enh_stop"))

common_pair <- subset(common_pair,common_pair$act==1&common_pair$dist>2e4)


gene_annotation <- read.table("../gene_annotation/gene_annotation_V19.txt")

housekeeping_gene <- read.table("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/HK_housekeeping_genes.txt")
housekeeping_gene <- as.character(housekeeping_gene[,1])

housekeeping_gene <- as.character(subset(gene_annotation[,7],gene_annotation[,6]%in%housekeeping_gene))

hk_gene <- read.csv("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/journal.pone.0022859.s008.csv",skip=1,header=T)[,2]

hk_gene <- as.character(hk_gene)

hk_gene <- as.character(subset(gene_annotation[,7],gene_annotation[,6]%in%hk_gene))

housekeeping_gene <- union(hk_gene,housekeeping_gene)


common_pair$hk_gene <- common_pair$gene_name%in%housekeeping_gene

common_pair$hk_rm_CaptureC <- common_pair$CaptureC * (common_pair$hk_gene==0)

common_pair <- na.omit(common_pair)

library(PRROC)

jeme.roc <- roc.curve(scores.class0=common_pair$jeme,weights.class0=common_pair$CaptureC,curve=T)
my.roc <- roc.curve(scores.class0=common_pair$pp,weights.class0=common_pair$CaptureC,curve=T)

my.pr <- pr.curve(scores.class0=common_pair$pp,weights.class0=common_pair$CaptureC,curve=T)

jeme.pr <- pr.curve(scores.class0=common_pair$jeme,weights.class0=common_pair$CaptureC,curve=T)

plot_roc <- function(x,main,col){
  
  plot(x$curve[,1],x$curve[,2],type="l",xlab="FPR",ylab="TPR",main=main,col=col,xlim=c(0,1),ylim=c(0,1))
  par(new=T)
  
}


pdf(paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/figures/systematic_conparison_updated/3_way_comparison/",ct,"_",gs_name,".pdf"))

plot_roc(jeme.roc,"","blue")
plot_roc(my.roc,"","red")


legend("bottomright",legend=c(paste0("Integrative model,AUC:",round(my.roc$auc,2)),paste0("JEME,AUC:",round(jeme.roc$auc,2))),col=c("red","blue"),lty=1)

par(new=F)

plot_pr <- function(x,main,col){
  
  plot(x$curve[,c(1,2)],type="l",xlim=c(0,1),ylim=c(0,y_lim),col=col,xlab = "Recall",ylab="Precision")
  
  par(new=T)
  
}

y_lim <- max(c(jeme.pr$curve[,2],my.pr$curve[,2]))

plot_pr(jeme.pr,"","blue")

plot_pr(my.pr,"","red")

legend("topright",legend=c(paste0("Integrative model,AUPR:",round(my.pr$auc.integral,5)),paste0("JEME,AUPR:",round(jeme.pr$auc.integral,5))),col=c("red","blue"),lty=1)

dev.off()

