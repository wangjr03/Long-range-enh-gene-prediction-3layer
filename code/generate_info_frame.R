args <- commandArgs(T)

ct <- as.numeric(args[1])



load("../data/gene_enh_frame.Rdata")

n <- dir(paste0("../output_exons_updated/",ct,"/i30_D5_DD10_C1_G200_K0.07_T10_B1_qt_pearson_NA/"))

n <- length(n)-1

load(paste0("../output_exons_updated/",ct,"/i30_D5_DD10_C1_G200_K0.07_T10_B1_qt_pearson_NA/iteration_",n,".Rdata"))

p_link <- sapply(post_p,function(x){sapply(x,function(y){
  y[1,]
})})

p_link.vec <- unlist(p_link)

load(paste0("../input_data_updated/",ct,"/input_data.Rdata"))


#prepare act information
pair_act <- c()
pot_act_long_list <- c()

for(i in 1:length(pot_act_list)){
  
  pot_act_long_list <- c(pot_act_long_list,pot_act_list[[i]])
  
}

pot_act_long_list <- sapply(pot_act_long_list,floor)

# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/No_conservation_motif/expressed_tf/input_data/gene_act")

gene_act <- unlist(gene_act)
gene_act <- floor(gene_act)
pair_act <- sapply(as.list(1:length(gene_act)),function(i){as.numeric(pot_act_long_list[[i]]*gene_act[i])>0})
pair_act <- unlist(pair_act)



# setwd(args[2])
# 
# ind <- list()
# for(i in 1:500){
#   
#   ind[[i]] <- data.table::fread(paste0(i,".txt"))
#   
# }
# 
# ind <- do.call(rbind,ind)
# capturec <- as.data.frame(ind)
# 
# 
# setwd("../")
# 
# 
# setwd(args[3])
# 
# ind <- list()
# for(i in 1:500){
#   
#   ind[[i]] <- data.table::fread(paste0(i,".txt"))
#   
# }
# ind <- do.call(rbind,ind)
# gtex <- as.data.frame(ind)


# setwd("../")


info_frame <- gene_enh_frame

info_frame$dist <- unlist(pot_dist_list)
info_frame$corr <- unlist(pot_corr_list)
info_frame$act <- pair_act
#info_frame$val <- info_frame$ind*info_frame$act
info_frame$tf_count <- unlist(tf_count_list)
info_frame$pp <- p_link.vec

# info_frame$CaptureC <- CaptureC[,1]
# 
# info_frame$GTEX <- gtex[,1]


#correct p with enh act
act_ind <- which(info_frame$act==1)


load("../data/enh_act_mat.Rdata")

info_frame$enh_act <- enh_act_mat[info_frame$enh,ct]


sub_info_frame <- subset(info_frame,info_frame$act==1)





print("#Calculating...")

L <- info_frame$pp

#d1 <- density(act[t_enh,cell_type],bw=0.5)
#d2 <- density(act[f_enh,cell_type],bw=0.5)
load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/non_linked_enhancer")
load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/linked_enhancer")



get_prob <- function(x,d){
  
  id <- which(d$x>x)[1]
  
  return(d$y[id])
  
  
}

get_prob <- Vectorize(get_prob,"x")

uniq_act <- sort(unique(info_frame$enh_act))

act_bf <- get_prob(uniq_act,d2)/get_prob(uniq_act,d1)

dic <- data.frame(uniq_act,act_bf)

id <- match(sub_info_frame$enh_act,dic[,1])

info_frame$act_bf <- 0

info_frame$act_bf[act_ind] <- dic[id,2]

info_frame$corrected_p <- 1/(1+(1/info_frame$pp-1)*info_frame$act_bf)

#load("..//")

dir.create(paste0("../final_prediction/",ct))

setwd(paste0("../final_prediction/",ct))

if(!"info_frame.Rdata"%in%dir()){
  
  save(info_frame,file="info_frame.Rdata")
  
}

