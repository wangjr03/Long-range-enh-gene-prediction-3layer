library(data.table)

args <- commandArgs(T)

index <- as.numeric(args[2])
enh_coord_path <- args[1]

file <- paste0("/mnt/ls15/scratch/users/wangha73/ENCODE_ROADMAP_127/RNA_Seq_bedGraph/",dir("/mnt/ls15/scratch/users/wangha73/ENCODE_ROADMAP_127/RNA_Seq_bedGraph/")[index])
names_order <- paste0("chr",c(1:22,"X"))

histone <- fread(file)
histone <- as.data.frame(histone)
histone[,1] <- factor(histone[,1],levels = names_order)
histone <- histone[order(histone[,1],histone[,2]),]

enh_coord <- read.table(enh_coord_path)
enh_coord[,2] <- factor(enh_coord[,2],levels = names_order)

enh_coord[,1] <- as.character(enh_coord[,1])

enh_coord <- enh_coord[,c(2,3,4)]

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

if(!"../gene_expr_promoter"%in%dir()){
  dir.create("../gene_expr_promoter")
}

save(average_score,score_chr,file=paste0("../gene_expr_promoter/",index,"_RNA_seq.Rdata"))
