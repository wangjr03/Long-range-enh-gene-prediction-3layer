#load data
library(PRROC)

args <- commandArgs(T)

cell_type <- args[1]

gs <- args[2]

gs_id <- as.numeric(args[3])

gs_name <- args[4]



setwd(paste0("../output_exons_updated/",cell_type))



load("info_frame.Rdata")

#load shuffle score
load("i30_D5_DD10_C1_G200_K0.07_T10_B1_qt_pearson_TF/iteration_7.Rdata")


p_link <- sapply(post_p,function(x){sapply(x,function(y){
  y[1,]
})})

p_link.vec <- unlist(p_link)

info_frame$shuffle_TF <- p_link.vec



#load 1 g score
load("i30_D5_DD10_C1_G2_K0.07_T10_B1_qt_pearson_NA/iteration_5.Rdata")


p_link <- sapply(post_p,function(x){sapply(x,function(y){
  y[1,]
})})

p_link.vec <- unlist(p_link)

info_frame$one_group <- p_link.vec



link_frame <- read.table(gs)

dist <- abs(link_frame[,4]-link_frame[,2])

dist <- dist[which(dist<2e6)]

d_breaks <- seq(0,2e6,length.out=200)
# bin <- apply(info_frame[1:10,],1,function(x){
#   
#   loc <- which(d_breaks>as.numeric(x[3]))[1]
#   if(loc>1&!is.na(loc)>0){
#     return(loc)}else{
#       return(0)
#     }
# }
# )

bin <- 1+ceiling(info_frame[,3]/10050.25)


count <- function(x){
  loc <- which(d_breaks>as.numeric(x))[1]
  if(loc>1&!is.na(loc)>0){
    return(loc)}else{
      return(0)
    }
}

count<- Vectorize(count)

link_frame_count <- count(dist)

bin_count <- as.numeric(table(factor(link_frame_count,levels=1:200)))

act_ind <- which(info_frame$act==1)


draw_pairs <- function(n){
  id <- c()
  for(i in 1:200){
    
    count <- floor(bin_count[i]*n/sum(bin_count))
    
    candidate <- which(bin==(i))
    
    id <- c(id,sample(intersect(candidate,act_ind),count))
    
  }
  
  return(id)
  
}

n_max = table(factor(bin[act_ind],levels=1:200))/bin_count*sum(bin_count)

n_max <- floor(min(n_max,na.rm=T))


draw_index <- draw_pairs(n_max)

roc.pp.list <- roc.shuffle.list <- roc.og.list <- c()

pr.pp.list <- pr.shuffle.list <- pr.og.list <- c()


for(i in 1:1000){
  print(i)
  draw_index_t <- sample(draw_index,1e4)
  
  roc.shuffle = roc.curve(scores.class0 = info_frame$shuffle_TF[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t])
  
  roc.og = roc.curve(scores.class0 = info_frame$one_group[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t])
  
  roc.pp = roc.curve(scores.class0 = info_frame$pp[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t])
  
  pr.shuffle = pr.curve(scores.class0 = info_frame$shuffle_TF[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t]>0)
  
  pr.og = pr.curve(scores.class0 = info_frame$one_group[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t]>0)
  
  pr.pp = pr.curve(scores.class0 = info_frame$pp[draw_index_t],weights.class0 = info_frame[,gs_id][draw_index_t]>0)
  
  roc.pp.list[i] <- roc.pp$auc
  
  roc.shuffle.list[i] <- roc.shuffle$auc
  
  pr.shuffle.list[i] <- pr.shuffle$auc.integral
  
  pr.pp.list[i] <- pr.pp$auc.integral
  
  roc.og.list[i] <- roc.og$auc
  
  pr.og.list[i] <- pr.og$auc.integral
  
  
}

library("ggpubr")


pdf(paste0("../../figures/systematic_conparison_updated/internal/ROC_",cell_type,"_",gs_name,".pdf"),width=5)

names=c("rand TF","one group","full model")

df <- data.frame(AUROC=c(roc.shuffle.list,roc.og.list,roc.pp.list),group=rep(names,each=1000))

df[,2] <- factor(as.character(df[,2]),levels = names)

#col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")


col=c("#4473A8","#70588E","#DB843D")

my_comparisons <- list(c(names[1],names[3]),c(names[2],names[3]))
ggboxplot(df, x = "group", y = "AUROC",col="group",palette=col)+ 
  stat_compare_means(comparisons = my_comparisons)

dev.off()


pdf(paste0("../../figures/systematic_conparison_updated/internal/PR_",cell_type,"_",gs_name,".pdf"),width = 5)

df <- data.frame(AUPR=c(pr.shuffle.list,pr.og.list,pr.pp.list),group=rep(names,each=1000))

df[,2] <- factor(as.character(df[,2]),levels = names)

col=c("#4473A8","#70588E","#DB843D")

my_comparisons <- list(c(names[1],names[3]),c(names[2],names[3]))

ggboxplot(df, x = "group", y = "AUPR",col="group",palette=col)+ 
  stat_compare_means(comparisons = my_comparisons)


dev.off()











