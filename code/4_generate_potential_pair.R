#2 parameters
#1: enhancer coordinates
#2: promoter coordinates

args <- commandArgs(T)

enh_coord <- read.table(args[1])

promoter <- read.table(args[2])

dist <- 2e6

enh_coord_c <- data.frame(chr=enh_coord[,1],center = floor((enh_coord[,2]+enh_coord[,3])/2))

promoter_c <- data.frame(gene_name=promoter[,1],chr=promoter[,2],center = floor((promoter[,3]+promoter[,4])/2))

promoter_c[,2] <- as.character(promoter_c[,2])

promoter_c[which(promoter_c[,2]=="chr23"),2] <- "chrX"

names_level <- paste0("chr",c(1:22,"X"))

enh_coord_c$chr <- factor(enh_coord_c$chr,levels = names_level)

promoter_c$chr <- factor(promoter_c$chr,levels = names_level)

enh_coord_c <- split(enh_coord_c,enh_coord_c$chr)

promoter_c <- split(promoter_c,promoter_c$chr)


pot_enh_list <- list()

for(i in 1:23){
  
  pot_enh_list[[i]] <- list()
  
  tmp_enh <- enh_coord_c[[i]]
  
  tmp_promoter <- promoter_c[[i]]
  
  n <- dim(tmp_promoter)[1]
  
  for(j in 1:n){
    
    tmp_gene <- tmp_promoter[j,]
    
    tmp_gene_name <- tmp_gene[,1]
    
    tmp_chr <- as.character(tmp_gene[,2])
    
    tmp_center <- tmp_gene[,3]
    
    pot_enh_list[[i]][[j]] <- which(tmp_enh$chr==tmp_chr&abs(tmp_enh$center-tmp_center)<dist)
    
  }
  
  names(pot_enh_list[[i]]) <- tmp_promoter[,1]
  
}

dir.create("../data")

save(pot_enh_list,file="../data/pot_enh_list.Rdata")
