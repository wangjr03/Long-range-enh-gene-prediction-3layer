#generate enh-promoter frame
load("../data/pot_enh_list.Rdata")
gene_name <- read.table("../promoter_activity/gene_name.txt",colClasses = "character")[,1]
enh_coord <- read.table("../consensus_enhancer_coord.txt")

names_level <- paste0("chr",c(1:22,"X"))

enh_coord[,1] <- factor(enh_coord[,1],levels = names_level)

enh_coord_count <- as.numeric(table(enh_coord[,1]))

enh_coord_count <- cumsum(c(0,enh_coord_count))


gene_enh_frame <- list()

k <- 1

for(i in 1:23){
  
  n <- length(pot_enh_list[[i]])
  
  for(j in 1:n){
    
    if(length(pot_enh_list[[i]][[j]])==0){next}
    
    tmp_name <- names(pot_enh_list[[i]])[j]
    
    tmp_id <- match(tmp_name,gene_name)
    
    tmp_enh_id <- enh_coord_count[i]+pot_enh_list[[i]][[j]]
    
    tmp_frame <- data.frame(gene=tmp_id,enh=tmp_enh_id)
    
    gene_enh_frame[[k]] <- tmp_frame
    
    k <- k+1
    
  }
  
}

gene_enh_frame <- do.call(rbind,gene_enh_frame)

save(gene_enh_frame,file="../data/gene_enh_frame.Rdata")