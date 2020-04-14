#get final predictions based on score cutoff or distance control
args <- commandArgs(T)

ct <- as.numeric(args[1])

load(paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/output_exons_updated/score/",ct,"_scores.Rdata"))

mode <- args[2]

enh_frame <- read.table("../consensus_enhancer_coord.txt")

promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")


#if apply distance control

if(mode=="score"){
  
  score <- as.numeric(args[3])
  
  dir.create(paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/output_exons_updated/final_prediction/",score))
  
  pred <- subset(info_frame,info_frame$corrected_p>score)
  
  pred_frame <- cbind(enh_frame[pred$enh,],promoter_frame[pred$gene,],pred[,c(3,4,10)])
  
  write.table(pred_frame,paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/output_exons_updated/final_prediction/",score,"/",ct,".txt"),
              col.names = F,row.names = F,sep="\t",quote = F)
  
  
}else if(mode=="distance"){
  
  score <- as.numeric(args[3])
  
  sub_info_frame <- subset(info_frame,info_frame$corrected_p>score)
  
  dir.create(paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/output_exons_updated/final_prediction/distance","_",score))
  
  dist <- read.table(args[4])
  
  sp_count <- as.numeric(args[5])
  
  dist <- abs(dist[,4]-dist[,2])
  
  dist <- subset(dist,dist < 2e6)
  
  d_breaks <- seq(0,2e6,length.out=201)
  
  count <- function(x){
    loc <- which(d_breaks>as.numeric(x))[1]
    if(loc>1&!is.na(loc)>0){
      return(loc)}else{
        return(0)
      }
  }
  
  count<- Vectorize(count)
  
  link_frame_count <- count(dist)
  
  bin_count <- as.numeric(table(factor(link_frame_count,levels=1:201)))
  
  
  draw_pairs <- function(bin_sample_count,bin){
    id <- list()
    for(i in 1:201){
      
      count <- ceiling(bin_sample_count[i])
      
      candidate <- which(bin==(i))
      
      candidate <- sub_info_frame[candidate,]
      
      if(length(candidate)==0|count==0) next
      
      candidate <- candidate[order(candidate[,10],decreasing=T),]
      
      id[[i]] <- candidate[1:count,]
      
    }
    
    id <- do.call(rbind,id)
    
    return(id)
    
  }
  
  my_bin <- ceiling(sub_info_frame[,3]/1e4)+1

    
    
  n_max = table(factor(my_bin,levels=1:201))/bin_count*sum(bin_count)
    
  n_max <- floor(min(n_max,na.rm=T))
  
  bin_sample_count <- bin_count/sum(bin_count)*sp_count  
  
  if(any(bin_sample_count>n_max)){
    
    print("Sampled too much! Use max sample in some bins")
    
    bin_sample_count[which(bin_sample_count>table(factor(my_bin,levels=1:201)))] <- table(factor(my_bin,levels=1:201))[which(bin_sample_count>table(factor(my_bin,levels=1:201)))] 
    
  }
    
  my_draw_index <- na.omit(draw_pairs(bin_sample_count,my_bin))
    
    
  pred_frame <- cbind(enh_frame[my_draw_index$enh,],promoter_frame[my_draw_index$gene,],my_draw_index[,c(3,4,10)])
  
  write.table(pred_frame,paste0("/mnt/gs18/scratch/users/wangha73/ENCODE_ROADMAP_127/output_exons_updated/final_prediction/distance","_",score,"/",ct,".txt"),
              col.names = F,row.names = F,sep="\t",quote = F)
  
}
