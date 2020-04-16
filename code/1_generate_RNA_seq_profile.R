#get score for enhancers
library(data.table)

args <- commandArgs(T)

index <- as.numeric(args[2])
enh_coord_path <- args[1]

RNA_seq_file <- args[3]


file <- paste0(args[3],dir(args[3])[index])
names_order <- paste0("chr",c(1:22,"X"))

histone <- fread(file)
histone <- as.data.frame(histone)
histone[,1] <- factor(histone[,1],levels = names_order)
histone <- histone[order(histone[,1],histone[,2]),]

enh_coord <- read.table(enh_coord_path,header=T)
enh_coord[,1] <- factor(enh_coord[,1],levels = names_order)

enh_coord[,4] <- as.character(enh_coord[,4])

#merge exons
# 
# merged_exons <- list()
# 
# tmp_gene_id <- enh_coord[1,7]
# 
# tmp_gene_name <- enh_coord[1,8]
# 
# tmp_start <- enh_coord[1,4]
# 
# tmp_stop <- enh_coord[1,5]
# 
# tmp_chr <- enh_coord[1,1]
# 
# k <- 1
# 
# for(i in 2:dim(enh_coord)[1]){
#   
#   if(enh_coord[i,7]!=tmp_gene_id){
#     
#     merged_exons[[k]] <- data.frame(tmp_chr,tmp_start,tmp_stop,tmp_gene_id,tmp_gene_name)
#     
#     tmp_gene_id <- enh_coord[i,7]
#     
#     tmp_gene_name <- enh_coord[i,8]
#     
#     tmp_start <- enh_coord[i,4]
#     
#     tmp_stop <- enh_coord[i,5]
#     
#     tmp_chr <- enh_coord[i,1]
#     
#     k <- k+1
#     
#   }else{ 
#     
#     if(enh_coord[i,4]>tmp_stop){
#     
#     merged_exons[[k]] <- data.frame(tmp_chr,tmp_start,tmp_stop,tmp_gene_id,tmp_gene_name)
#     
#     tmp_gene_id <- enh_coord[i,7]
#     
#     tmp_gene_name <- enh_coord[i,8]
#     
#     tmp_start <- enh_coord[i,4]
#     
#     tmp_stop <- enh_coord[i,5]
#     
#     tmp_chr <- enh_coord[i,1]
#     
#     k <- k+1
#     
#   }else{
#     
#     if(enh_coord[i,5]>tmp_stop){
#     
#     tmp_stop <- enh_coord[i,5]
#     
#     }
#   }
#   }
#   
# }
# 
# merged_coord <- do.call(rbind,merged_exons)

enh_coord <- enh_coord[,c(1,2,3)]

#split data by chromosome
enh_coord_split <- split(enh_coord,enh_coord[,1])

histone <- split(histone,histone[,1])

score_chr <- list()

average_score <- list()

for(chr in 1:23){
  
  temp_histone <- histone[[chr]]
  
  temp_enh_coord <- enh_coord_split[[chr]]
  
  j <- 1
  
  score <- list()
  
  score[[1]] <- logical(0)
  
  for(i in 1:dim(temp_histone)[1]){
    
    if(temp_histone[i,1]!=temp_enh_coord[j,1]){break}
    
    if(!temp_histone[i,3]>temp_enh_coord[j,2]){
      
      next
      
    }else if(!temp_histone[i,2]<temp_enh_coord[j,3]){
      
      j <- j+1
      print(j)
      if(j > dim(temp_enh_coord)[1]){break}
      score[[j]] <- logical(0)
      
    }
    
    
    if(!temp_histone[i,2]>temp_enh_coord[j,3]){
    score[[j]] <- rbind(score[[j]],data.frame(id=i,score=temp_histone[i,4]))}
    
    
    
  }
  average_score[[chr]] <- sapply(score, function(x){
    if(length(x)>0){
      mean(x[,2])
    }else{
      
      0
      
    }
    
    }
    )
  score_chr[[chr]] <- score
  
}

average_score <- unlist(average_score)

if(!"../gene_expr_exons"%in%dir()){
  dir.create("../gene_expr_exons")
}

save(average_score,score_chr,file=paste0("../gene_expr_exons/",index,"_RNA_seq.Rdata"))
