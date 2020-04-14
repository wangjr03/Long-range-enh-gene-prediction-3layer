args <-commandArgs(T)

cell_type <- as.numeric(args[1])

enh_coord <- read.table("../consensus_enhancer_coord.txt")
load("../data/enh_act_mat.Rdata")

enh_id <- enh_coord

enh_act <- enh_act_mat

enh_loc <- prodlim::row.match(enh_coord,enh_id)

enh_act <- enh_act[enh_loc,cell_type]

enh_center <- (enh_coord[,2]+enh_coord[,3])/2

enh_act <- data.frame(chr=enh_coord[,1],center=floor(enh_center),act=enh_act)

enh_act$act <- floor(enh_act$act)

#gene_act
gene_name <- read.table("../promoter_activity/gene_name.txt",colClasses = "character")[,1]

gene_name <- sapply(strsplit(gene_name,"[.]"),function(x){x[[1]][1]})

load("../data/gene_act_mat_promoter.Rdata")

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ref <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"), filters = "ensembl_gene_id", values = gene_name,mart = mart)

trs_loc <- match(ref[,1],gene_name)

gene_act <- data.frame(trascript_id=ref[,2],act=gene_act_mat[trs_loc,cell_type])

gene_act[,2] <- floor(gene_act[,2])

gene_act <- subset(gene_act,gene_act[,2]>0)

enh_coord <- subset(enh_coord,enh_act$act>0)

enh_act <- subset(enh_act,enh_act$act>0)

pth <- paste0("../IM-PET/",cell_type,"/")

dir.create(pth,recursive = T)

setwd(pth)

write.table(enh_coord,"enh_coord",col.names=F,row.names = F,sep="\t",quote=F)
write.table(enh_act,"enh_act",col.names=F,row.names = F,sep="\t",quote=F)
write.table(gene_act,"gene_act",col.names=F,row.names = F,sep="\t",quote=F)
