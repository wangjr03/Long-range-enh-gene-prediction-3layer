#load enhancer/promoter activity profile
calculate_dist <- function(x,y){
  
  promoter_c <- floor((tmp_gene[x,3]+tmp_gene[x,4])/2)
  
  enh_c <- floor((tmp_enh[y,2]+tmp_enh[y,3])/2)
  
  return(abs(promoter_c-enh_c))
  
}

calculate_dist <- Vectorize(calculate_dist,vectorize.args = "y")

enh_coord <- read.table("../consensus_enhancer_coord.txt")

promoter <- read.table("../promoter_activity/promoter_frame.txt")

promoter <- read.table("../promoter_activity/promoter_frame.txt")

promoter[,2] <- as.character(promoter[,2])

promoter[which(promoter[,2]=="chr23"),2] <- "chrX"

names_level <- paste0("chr",c(1:22,"X"))

enh_coord[,1] <- factor(enh_coord[,1],levels = names_level)

promoter[,2] <- factor(promoter[,2],levels = names_level)

#get act profile in each chromosome
enh_coord_chr <- split(enh_coord,enh_coord[,1])

promoter_chr <- split(promoter,promoter[,2])
#load all pot enhancer promoter pool
load("../data/pot_enh_list.Rdata")

pot_dist_list <- list()

for(i in 1:23){
  
  tmp_enh <- enh_coord_chr[[i]]
  
  tmp_gene <- promoter_chr[[i]]
  
  tmp_pot_enh <- pot_enh_list[[i]]
  
  n <- length(tmp_pot_enh)
  
  pot_dist_list[[i]] <- sapply(1:n,function(x){
    
    calculate_dist(x,tmp_pot_enh[[x]])
    
  })
  
  
}

save(pot_corr_list,file="../data/pot_dist_list.Rdata")



