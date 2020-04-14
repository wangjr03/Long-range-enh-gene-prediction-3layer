all_start_time <- Sys.time()

#args <- c("5","10","1","200","test","0.07",
#          "/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/input_data/input_data/input_data.Rdata",
#          "1","ALL")

args <- commandArgs(TRUE)
D <- as.numeric(args[1])  ## KL coefficient for enh-gene interaction, initial number
DD <- as.numeric(args[2])  ## KL coefficient for gene group
C <- as.numeric(args[3])  ## correlation coefficient
B <- 1 ## distance coefficient 
D_range <- c(0,30) #range of D, should from seed
Temperature <- 10
d_m <- as.numeric(args[8])
input_data <- args[7]
############################## simulated annealing parameter #######################
K <- as.numeric(args[6]) #tolerence of correlation kl divergence
#T should in the same order with e
###################################################################################

randomize_term <- args[9]
num_of_group <- as.numeric(args[4])   ## how many gene group (G)
dir_path <- args[5]
post_p <- list()
#cell <- 13    ## which cell type  "CD4_Memory_Primary_Cells"
num_iter <- 100  ## how many iteration
num_cores <- 28  ## how many cores used for computing
#num_of_group <- 500   ## how many gene group
pcount_gene <- 3 
pcount_corr <- 1
print(paste("num_iter=", num_iter, sep=""))
print(paste("num_of_group=", num_of_group, sep=""))
print(paste("D=", D, sep = ""))
print(paste("DD=", DD, sep = ""))
print(paste("C=", C, sep = ""))
print(paste("K=", K, sep = ""))
print(paste("Randomize=", randomize_term, sep = ""))
print(paste("d_m=", d_m, sep = ""))



library(parallel)
library(entropy)
library(pheatmap)
library(LaplacesDemon)
library(truncnorm)
## get chromosome name
library(data.table)
enh <- fread("../consensus_enhancer_coord.txt")
enh <- as.data.frame(enh)
chrname <- unique(enh$V1)
length <- length(chrname) ## how many chromosome
chr_num <- c(1:length)

## load the purified data
# load("/mnt/research/compbio/wanglab/hongjie/all_TF_GE.Rdata")
# load("/mnt/research/compbio/wanglab/hongjie/final_input_data.Rdata")
# load("/mnt/research/compbio/wanglab/hongjie/initial_group_norm.Rdata")
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/corr_init_final_input.Rdata")
# load("/mnt/research/compbio/wanglab/hongjie/pot_distance_list.Rdata")
# 
# #load quantile normalized correlation
# load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/pot_corr_list_10ct_spearman_1.8M")
# pot_corr_list <- pot_corr_list_10ct_spearman


load(input_data)
pot_distance_list <- pot_dist_list





#randomize a list, by keeping the number in each element
randomize_list <- function(x){
  
  n <- length(x)
  
  count <- sapply(x,length)
  
  rand_ele <- sample(unlist(x))
  
  rand_list <- list()
  
  for(i in 1:n){
    
    rand_id <- sample(1:length(rand_ele),count[i])
    
    rand_list[[i]] <- rand_ele[rand_id]
    
    rand_ele <- rand_ele[-rand_id]
    
    
  }
  
  return(rand_list)
  
}

################# select which one to randomize #########################################
if(randomize_term == "TF"){
  
  enh_tf_list <- sapply(enh_tf_list,function(x){
    
    apply(x,2,sample)
    
  })
  
}else if(randomize_term == "corr"){
  
  pot_corr_list <- sapply(pot_corr_list,randomize_list)
  
}else if(randomize_term == "dist"){
  
  pot_distance_list <- sapply(pot_distance_list,randomize_list)
  
}


################################shuffle TF on all enhancers#############################
if(randomize_term=="ALL"|randomize_term=="bg"){
  
  names_order <- paste0("chr",c(1:22,"X"))
  
  print("shuffling TF matrix across all enhancers")
  #generate a enh_tf matrix for all enhancers
  #setwd("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/hongjie/data/enhancer_split/")
  
  all_enh <- list()
  all_enh_tf <- list()
  
  for(i in 1:23){
    if(i != 23){
      
      all_enh[[i]] <- subset(enh,enh[,1]==names_order[i])
      
      all_enh_tf[[i]] <- matrix(0,nrow=dim(all_enh[[i]])[1],ncol=dim(enh_tf_list[[i]])[2])
      
      all_enh_tf[[i]][use_pos_list[[i]],] <- enh_tf_list[[i]]
      
    }else{
      
      all_enh[[i]] <- subset(enh,enh[,1]==names_order[i])
      
      all_enh_tf[[i]] <- matrix(0,nrow=dim(all_enh[[i]])[1],ncol=dim(enh_tf_list[[i]])[2])
      
      all_enh_tf[[i]][use_pos_list[[i]],] <- enh_tf_list[[i]]
      
    }
    
    
  }
  
  #shuffle the all_enh_tf matrix
  
  
  all_enh_tf <- sapply(all_enh_tf,function(x){
    
    apply(x,2,sample)
    
  })
  
  
  #generate new use_pos
  
  use_pos_list <- sapply(all_enh_tf,function(x){
    
    count <- apply(x,1,sum)
    
    return(which(count>0))
    
  })
  
  #get count_list 
  count_list <- sapply(all_enh_tf,function(x){
    
    count <- apply(x,1,sum)
    
    return(count)
    
  })
  
  #get tf_count_list
  
  
  for(i in 1:23){
    
    tf_count_list[[i]] <- lapply( pot_enh_list[[i]],function(x){
      
      return(count_list[[i]][x])
      
    })
    
  }
  
  #extract enh_tf_matrix
  
  enh_tf_list <- sapply(all_enh_tf,function(x){
    
    id <- which(apply(x,1,sum)>0)
    
    return(x[id,])
    
  })
  
  if(randomize_term=="bg"){
    
    print("shuffing correlation")
    pot_corr_list <- sapply(pot_corr_list,randomize_list)
    #print("shuffling distance")
    #pot_distance_list <- sapply(pot_distance_list,randomize_list)
  }
  
  setwd(dir_path)
  save(tf_count_list,file="shuffeled_tf_list")
  save(use_pos_list,file="shuffeled_use_pos_list")
  save(enh_tf_list,file="shuffeled_enh_tf_list")
  save(pot_corr_list,file="shuffeled_pot_corr_list")
  #save(pot_distance_list,file="shuffeled_pot_distance_list")
  
  
  
  
}

group_infor <- as.vector(cutree(tree_gene, k = num_of_group))

# file_name_iter <- paste("/mnt/research/compbio/wanglab/hongjie/output/iteration_result_1115/D1DD1C1_KL/iteration_",1:num_iter, ".Rdata",sep="")
file_name_iter <- paste(dir_path,"/iteration_",1:num_iter, ".Rdata",sep="")

all_gene_name_list <- list()
for (g in 1:nrow(gene_infor_matrix)) {
  all_gene_name_list[[g]] <- as.character(gene_infor_matrix[g,1])
}

################### randomly assign enhancer-gene interaction ######################
# random_assign_enh <- function (pot) {
#   if (length(pot) >= 10) {
#     pos <- sample(1:length(pot), size = 10)
#     inter <- pot[pos]
#   } else {
#     inter <- pot
#   }
#   return(inter)
# }
# 
# ## update inter_enh_list
# random_inter_enh_list <- list()
# for (i in 1:23) {
#   all_enh_pot <- pot_enh_list[[i]]
#   new_inter <- lapply(all_enh_pot, random_assign_enh)
#   random_inter_enh_list[[i]] <- new_inter
#   
#   if (class(new_inter) != "list") {
#     print(paste("chr",i,"not list"))
#   }
# }
# 
# ## update inter_corr_list 
# extract_corr <- function(pot_corr, inter, pot) {
#   pos <- which(pot %in% inter)
#   inter_corr <- pot_corr[pos]
#   return(inter_corr)
# }
# 
# new_inter_corr_list <- list()
# for (i in 1:23) {
#   pot_corr <- pot_corr_list[[i]]
#   inter <- inter_enh_list[[i]]
#   pot <- inter_enh_list[[i]]
#   
#   new_corr <- mapply(extract_corr, pot_corr, inter, pot)
#   new_inter_corr_list[[i]] <- new_corr
# }
# 
# ## update TF_GE_list
# new_TF_GE_list <- list()
# for (i in 1:23) {
#   gene_name_all_exp <- gene_name_list[[i]]
#   num_gene_exp <- length(gene_name_all_exp)
#   enh_tf <- enh_tf_list[[i]]
#   num_of_tf <- ncol(enh_tf)
#   all_enh_inter_exp <- inter_enh_list[[i]]
#   
#   TF_GE <- matrix(0,num_gene_exp,num_of_tf)
#   for (g in 1:num_gene_exp) {
#     sub_enh_pos <- all_enh_inter_exp[[g]]
#     sub_enh_tf <- matrix(enh_tf[sub_enh_pos,], length(sub_enh_pos), num_of_tf)
#     
#     TF_GE[g, ] <- as.vector(apply(sub_enh_tf, 2, mean))
#   }
#   #rownames(TF_GE) <- gene_name_all_exp
#   TF_GE_df <- data.frame(name = gene_name_all_exp, TF_GE)
#   colnames(TF_GE_df) <- c("name",colnames(enh_tf))
#   
#   new_TF_GE_list[[i]] <- TF_GE_df
# }
#   TF_GE_list <- new_TF_GE_list
# 
# ## update all_TF_GE 
# pre_all_TF_GE <- all_TF_GE
# #### construct Gene-TF matrix for all chromosome 
# replace_na <- function(vector) {
#   pos <- which(is.na(vector) == TRUE)
#   vector[pos] <- 0
#   return(vector)
# }
# 
# all_TF_GE <- TF_GE_list[[1]]
# for (i in 2:23) {
#   new <- TF_GE_list[[i]]
#   #new <- apply(new, 1, normalize)
#   all_TF_GE <- merge(all_TF_GE, new, all=TRUE)
# }
# 
# all_gene_list <- all_TF_GE[,1]
# all_TF_GE <- as.matrix(all_TF_GE[,-1])
# all_TF_GE <- apply(all_TF_GE, 2, as.numeric)
# all_TF_GE <- apply(all_TF_GE, 2, replace_na)
# row.names(all_TF_GE) <- all_gene_list
# 
# #### reorder the gene in all_TF_GE according gene_infor_matrix
# pos_gene_match <- match(gene_infor_matrix[,1], rownames(all_TF_GE))
# all_TF_GE <- all_TF_GE[pos_gene_match,]
# 
# #### reorder TF in every TF_GE matrix
# for (i in 1:23) {
#   TF_GE <- all_TF_GE[which(gene_infor_matrix[,2] == i),]
#   TF_GE_list[[i]] <- TF_GE
# }
# 
# ## update group membership initialization
# library(pheatmap)
# data <- all_TF_GE
# data <- data[,-as.vector(which(apply(data,2,sum) == 0))]
# phm <- pheatmap(data,color = colorRampPalette(c("snow", "red"))(50),clustering_distance_row = "correlation", clustering_distance_col = "correlation", fontsize=9, fontsize_row=1, fontsize_col=3)
# tree_gene <- phm$tree_row
# tree_tf <- phm$tree_col
# group_infor <- as.vector(cutree(tree_gene, k = 200))
# 
# ## save data
# save(
#   all_TF_GE,
#   file = "/mnt/research/compbio/wanglab/hongjie/random_all_TF_GE.Rdata"
# )
# 
# save(
#   all_TF_GE,
#   TF_GE_list,
#   enh_tf_list,
#   tf_name_list,
#   all_tf_name,
#   
#   use_var_num,
#   use_pos_list,
#   
#   gene_name_list,
#   gene_infor_matrix,
#   pot_enh_list,
#   inter_enh_list,
#   pot_corr_list,
#   inter_corr_list,
#   
#   group_infor,
#   tree_gene,
#   tree_tf,
#   
#   file = "/mnt/research/compbio/wanglab/hongjie/random_final_input.Rdata"
# )

######################## define basic functions ###########################

#linear decrease, AUC 0.93
# calculate_D <- function(x){
#   
#   if(x<4.5){return(1)}else if(x>5){
#     
#     return(0)
#     
#   }else{
#     
#     -2*x+10
#     
#   }
#   
#   
# }

#polynorm decrease

calculate_D <- function(x){
  
  if(x<4.5){return(1)}else if(x>5){
    
    return(1)
    
  }else{
    
    return(1)
    #-4*x^2+36*x-80
    
  }
  
  
}



## function: to calculate KL
KL <- function(p,q) {
  p <- p/(sum(p) + 10^(-60))
  q <- q/(sum(q) + 10^(-60))
  p <- p + 10^(-60)
  q <- q + 10^(-60)
  kl <- min(KL.empirical(p,q), KL.empirical(q,p))
  return(kl)
  #(1-15)return(KLD(p,q)[[7]])
}

## function: to calculate KL only on none-zero dimensons
KL_indic <- function(p,q) {
  p <- p/(sum(p) + 10^(-60))
  q <- q/(sum(q) + 10^(-60))
  q <- q + 10^(-60)
  indicator <- (q > 0)
  kl <- sum(KLD(p,q)[[2]] * indicator)
  return(kl)
}

## function: to split names in one ""
spl_name <- function(one_name) {
  aft_spl <- unlist(strsplit(as.character(one_name), ","))
  return(aft_spl)
  ## end
}

## function: to replace NA with 0 
replace_na <- function(vector) {
  pos <- which(is.na(vector) == TRUE)
  vector[pos] <- 0
  return(vector)
}

## function: to replace Inf with 1
replace_inf <- function(vector) {
  pos <- which(vector == Inf)
  vector[pos] <- 1
  return(vector)
}

## function: To remove zero value in a vector
remove_zero <- function(x) {
  pos <- which(x == 0)
  y <- x[-pos]
  if (length(y) == 0 ) {
    y <- NULL
  }
  return(y)
}

## function: To normalize and get z-score
z_score <- function(x) {
  miu <- mean(x)
  sigma <- sd(x)
  new_x <- (x - miu) / sigma
  return(new_x)
}

## function: To get the proportion for each element of a vector
prop <- function(x){
  all <- sum(x)
  new <- x / all
  return(new)
}



## function: to decide whether V1 is in V2, and return 1 for TRUE 0 for FALSE
## usually V2 is a subset of V1
in_or_not <- function(x,y) {
  result <- ifelse(x %in% y, 1, 0)
  return(result)
}

## function: to decompose bin_link_list (a list of lists of vectors)
decomp_bin <- function(l) {
  v <- sapply(l, unlist)
  v <- sapply(v, as.vector)
  v <- as.vector(unlist(v))
}


## function: calculate the mean after add new observation
add_one_mean <- function(old_mean, num, obs) {
  
  summ <- old_mean*num + obs
  new_mean <- summ/(num+1)
  return(new_mean)
  
}

## function: calculate the mean after remove new observation
remove_one_mean <- function(old_mean, num, obs) {
  if (num > 1) {
    summ <- old_mean*num - obs
    new_mean <- summ/(num-1)
    neg_pos <- which(new_mean < 0) 
    if (length(neg_pos) > 0) {
      new_mean[neg_pos] <- 0
    }
    return(new_mean)
  } else {
    return(old_mean - old_mean)
  }
}

## fuction: to transform a matrix into a list(each row/column becomes a element)
trans_to_list <- function(m, ind) {
  new_list <- list()
  
  if (ind == 1) {
    for (d in 1:dim(m)[1]) {
      vec <- m[d,]
      new_list[[d]] <- vec
    }
  } else {
    for (d in 1:dim(m)[2]) {
      vec <- m[,d]
      new_list[[d]] <- vec
    }
  }
  
  return(new_list)
}


## function: to take the mean of each element in a list
take_mean_in_list <- function(m) {
  one_vector <- apply(m,2,mean)
  return(one_vector)
}


## function: normalize a probabilities vector and make their sum equal to 1
normalize_to_one <- function(vector) {
  sum <- sum(vector)
  new_vector <- vector/sum
  
  summ <- sum(new_vector)
  diff <- summ - 1
  
  if (diff < 0) {
    length <- length(vector)
    ave_diff <- diff/length
    new_vector <- new_vector - ave_diff
  } else if (diff > 0) {
    
    summ <- sum(new_vector)
    while (summ > 1) {
      none_zero_pos <- which(new_vector != 0)
      none_zero_vec <- new_vector[none_zero_pos]
      
      min_pos <- which.min(none_zero_vec)[1]
      none_zero_vec[min_pos] <- 0
      
      new_vector[none_zero_pos] <- none_zero_vec 
      summ <- sum(new_vector)
    }
    
    length <- length(vector)
    ave_diff <- (1 - sum(new_vector)) / length
    new_vector <- new_vector + ave_diff
    
  }
  
  return(new_vector)
}


## function: To normalize log-probabilities matrix
## each row is a possible outcome; each coloum is a part of probability
## return a matrix
normalize_log_prob <- function(log_probs) {
  log_probs <- as.matrix(log_probs)
  
  probs_matrix <- NULL
  for (l in 1:dim(log_probs)[2]) {
    
    vec_log_probs <- log_probs[ ,l]
    
    vec_prob <- rep(0, length(vec_log_probs))
    for (e in 1:length(vec_log_probs)) {
      one_log_prob <- vec_log_probs[e]
      deno <- exp(vec_log_probs - one_log_prob) 
      deno <- sum(deno) 
      prob <- 1 / deno
      vec_prob[e] <- prob
    }
    
    probs_matrix <- cbind(probs_matrix, vec_prob)
    
  }
  
  return(probs_matrix)
  
}


## function: normalize the log-probabilities and do the sampling
norm_and_samp <- function(vec) {
  probs_sep <- as.vector(normalize_log_prob(vec))
  probs <- normalize_to_one(probs_sep)
  samp_result <- sample(c(1:length(vec)), size = 1, prob = probs)
  return(samp_result)
}

samp_one <- function (vec) {
  result <- sample(c(1:length(vec)), size = 1, prob = vec)
  return(result)
}

## function: merge all matrixs in a list together
## all matrixs need to have the same column number
merge_matrix_list <- function(ls) {
  all <- NULL
  for (nl in 1:length(ls)) {
    all <- rbind(all, ls[[nl]])
  }
  return(all)
}


## function: to estimate the distribution and calculate the occurance probability for a sample
log_prob_occur <- function(data, one_sample) {
  
  data <- as.matrix(data)
  if (nrow(data) < 2) {
    data <- rbind(data, matrix(0, 2, ncol(data)))
  }
  
  all_log_prob <- NULL
  for (d in 1:dim(data)[2]) {
    est_density <- density(data[,d])
    x_pos <- which(est_density$x >= min(max(est_density$x),as.numeric(one_sample[d])))[1]
    width <- est_density$x[2] - est_density$x[1] 
    log_prob <- log( width * (10^(-20)+est_density$y[x_pos]) )
    all_log_prob <- c(all_log_prob, log_prob)
  }
  
  return(sum(all_log_prob))
  
}


## function: to get the estimated distributions for each dimensions
get_log_dist <- function(data) {
  data <- as.matrix(data)
  if (nrow(data) < 2) {
    data <- rbind(data, matrix(0, 2, ncol(data)))
  }
  
  all_dist <- list()
  for (d in 1:dim(data)[2]) {
    est_density <- density(data[,d])
    all_dist[[d]] <- est_density
  }
  
  return(all_dist)
}


## function: to get  the occurance probability for a new sample from a estimated distribution
get_log_prob <- function(all_dist, one_sample) {
  
  all_log_prob <- NULL
  for (d in 1:length(all_dist)) {
    est_density <- all_dist[[d]]
    x_pos <- which(est_density$x >= min(max(est_density$x),as.numeric(one_sample[d])))[1]
    width <- est_density$bw
    #width <- est_density$x[2] - est_density$x[1] 
    #width <- 1
    log_prob <- log( width * (10^(-20)+est_density$y[x_pos]) )
    all_log_prob <- c(all_log_prob, log_prob)
  }
  return(sum(all_log_prob))
  
}

## function: calculate accept probability of simulated annealing
accept_prob <- function(in_e,in_e_star,Temperature){
  
  p <- min(1,exp((in_e-in_e_star)/Temperature))
  
  
  return(p)
}

## function: decide accept or reject

decision <- function(in_p){
  
  sample_p <- runif(1)
  
  d <- ifelse(sample_p<in_p,1,0)
  
  return(d)
  
}


#function: calculate energy
calculate_energy <- function(in_group_tf_mean_data,in_group_corr_data,in_group_non_corr_data,K){
  
  temp_e <- sum(sapply(in_group_tf_mean_data,function(x){entropy(as.vector(x))}))
  
  y_non_corr <- density(unlist(in_group_non_corr_data),from=-1,to=1,n=512)$y
  
  y_corr <- density(unlist(in_group_corr_data),from=-1,to=1,n=512)$y
  
  temp_kl <- max(1,exp(K-KL(y_non_corr,y_corr)))
  
  e <- temp_e*temp_kl
  
  return(e)
}

######################## assign initial value for gene group ###########################

# library(pheatmap)
# #num_of_group <- K
# 
# data <- all_TF_GE
# data <- data[,-as.vector(which(apply(data,2,sum) == 0))]
# phm <- pheatmap(data,color = colorRampPalette(c("snow", "red"))(50),clustering_distance_row = "correlation", clustering_distance_col = "correlation", fontsize=9, fontsize_row=1, fontsize_col=3)
# tree_gene <- phm$tree_row
# tree_tf <- phm$tree_col
# group_infor <- as.vector(cutree(tree_gene, k = num_of_group))
# 
# # corr_matrix <-cor(t(all_TF_GE))
# # di <- dist(corr_matrix)
# # hc <- hclust(di, method = "complete")
# # group_infor <- as.vector(cutree(hc, k = num_of_group))
# 
# save(
#   group_infor,
#   tree_gene,
#   tree_tf,
#   file = "/mnt/research/compbio/wanglab/hongjie/initial_group.Rdata"
# )


#########################################################################################



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##
## iteration starts here
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
tf_sum <- apply(all_TF_GE, 2, sum) #(11/21)
count_TF_GE <- all_TF_GE #(11/21)
all_TF_GE <- apply(all_TF_GE, 2, prop) 
all_TF_GE <- apply(all_TF_GE, 2, replace_na)

## construct list to record the changes of each enhancer-gene pair
pre_bin_link_list <- list()
bin_link_list <- list()
for (i in 1:length) {
  all_enh_pot <- pot_enh_list[[i]]
  all_enh_inter <- inter_enh_list[[i]]
  bin_link_list[[i]] <- mapply(FUN=in_or_not, x=all_enh_pot, y=all_enh_inter)
}

## extract correlation data for non-interacted enhancer-gene pair
non_inter_corr_list <- list()
for (i in 1:length) {
  one_chr <- list()
  for (j in 1:length(pot_enh_list[[i]])) {
    
    #0729 edit 
    inter_pos <- which(pot_enh_list[[i]][[j]] %in% inter_enh_list[[i]][[j]])
    act_pos <- which(pot_act_list[[i]][[j]]*gene_act[[i]][[j]]>0)
    tf_pos <- which(tf_count_list[[i]][[j]]>0)
    
    #0905 relax to active enhancers
    non_link_pos <- setdiff(act_pos,inter_pos)
    
    if (length(non_link_pos != 0)) {
      one_gene <- pot_corr_list[[i]][[j]][non_link_pos]
    } else {
      if(length(act_pos)==0){
        
        one_gene <- logical(0)
        
      }else{
        
        #if all enhancers are linked, null correlation
        one_gene <- logical(0)
        
      }
      
    }
    one_chr[[j]] <- one_gene
  }
  non_inter_corr_list[[i]] <- one_chr
}

## construct distance distributions
inter_distance <- get_log_dist(all_inter_dist)
back_distance <- get_log_dist(unlist(pot_distance_list))


#######################################################initialize gene group data########################


## prepare gene group list
group_infor_list <- list()
for (g in 1:length(group_infor) ) {
  group_infor_list[[g]] <- group_infor[g]
}


## prepare TF data  for each gene group
mean_gene_tf <- apply(all_TF_GE, 2, mean)
group_tf_data <- list()
for (k in 1:num_of_group) {
  
  gene_pos <- which(group_infor == k)
  if (length(gene_pos) > 0) {
    
    if (length(gene_pos) > 1) {
      tf_for_this_group <- all_TF_GE[gene_pos, ]
    } else {
      tf_for_this_group <- rbind(all_TF_GE[gene_pos, ], mean_gene_tf)
    }
    group_tf_data[[k]] <- tf_for_this_group
    
  } else {
    #group_tf_data[[k]] <- matrix(0, pcount_gene, dim(all_TF_GE)[2])
    group_tf_data[[k]] <- matrix(rep(mean_gene_tf, pcount_gene), pcount_gene, dim(all_TF_GE)[2], byrow = T)
  }
}



group_tf_mean_data <- lapply(group_tf_data, take_mean_in_list)


# (11/15)
# ## prepare TF distribution for each gene group
# group_tf_dist <- list()
# for (k in 1:num_of_group) {
#   group_tf<- group_tf_data[[k]]
#   group_tf_dist[[k]] <- get_log_dist(group_tf)
# }


## prepare correlation data for each gene group
mean_gene_corr <- mean(unlist(inter_corr_list))
group_corr_data <- list()
for (k in 1:num_of_group) {
  gene_pos <- which(group_infor == k)
  
  if (length(gene_pos) != 0 ) {
    
    corr_one_group <- NULL
    for (l in gene_pos) {
      which_chr <- gene_infor_matrix[l,2]
      pos_in_chr <- gene_infor_matrix[l,3]
      corr_one_gene_inter <- inter_corr_list[[which_chr]][[pos_in_chr]]
      corr_one_group <- c(corr_one_group, corr_one_gene_inter)
    }
    if (length(corr_one_group) == 1) {
      corr_one_group <- c(corr_one_group, mean_gene_corr)
    }
    group_corr_data[[k]] <- as.matrix(corr_one_group)
    
  } else {
    #group_corr_data[[k]] <- rep(0, pcount_corr)
    group_corr_data[[k]] <- as.matrix(rep(mean_gene_corr, pcount_corr))
  }
  
}

## prepare non-intercted correlation data for each gene group
mean_non_gene_corr <- mean(unlist(non_inter_corr_list))
group_non_corr_data <- list()
for (k in 1:num_of_group) {
  gene_pos <- which(group_infor == k)
  
  if (length(gene_pos) != 0 ) {
    
    corr_one_group <- NULL
    for (l in gene_pos) {
      which_chr <- gene_infor_matrix[l,2]
      pos_in_chr <- gene_infor_matrix[l,3]
      corr_one_gene_non_inter <- non_inter_corr_list[[which_chr]][[pos_in_chr]]
      corr_one_group <- c(corr_one_group, corr_one_gene_non_inter)
    }
    if (length(corr_one_group) == 1) {
      corr_one_group <- c(corr_one_group, mean_non_gene_corr)
    }
    group_non_corr_data[[k]] <- as.matrix(corr_one_group)
    
  } else {
    #group_corr_data[[k]] <- rep(0, pcount_corr)
    group_non_corr_data[[k]] <- as.matrix(rep(mean_non_gene_corr, pcount_corr))
  }
  
}


## prepare correlation distribution for each gene group
group_corr_dist <- list()
data <- unlist(inter_corr_list)
d <- get_log_dist(data)
for (k in 1:num_of_group) {
  one_group_dist <- d
  group_corr_dist[[k]] <- one_group_dist
}

## prepare non-interacted correlation distribution for each gene group
group_non_corr_dist <- list()
data <- unlist(non_inter_corr_list)
d <- get_log_dist(data)
for (k in 1:num_of_group) {
  one_group_dist <- d
  group_non_corr_dist[[k]] <- one_group_dist
}







################################################### iteration starts ###############################################
for (n in 1:num_iter) {
  iter_start_time <- Sys.time()
  
  ######################## update enhancer-gene interaction ###########################
  
  enh_start <- Sys.time() ## update enhancer-gene interaction start time
  
  
  ################################################ remember previous result #################################
  pre_enh_list <- inter_enh_list
  pre_bin_link_list[[i]] <- bin_link_list[[i]]
  pre_post_p <- post_p
  mean_group_tf <- apply(all_TF_GE,2,mean)
  ## update interaction chr-by-chr
  for (i in 1:length) {
    
    # start_time <- Sys.time()
    print(i)
    ## Get data for this chromosome
    TF_GE <- TF_GE_list[[i]]
    all_tf_count <- tf_count_list[[i]]
    enh_tf <- enh_tf_list[[i]]  
    all_gene_name <- gene_name_list[[i]]
    all_enh_pot <- pot_enh_list[[i]]
    all_enh_inter <- inter_enh_list[[i]]
    all_dist_pot <- pot_distance_list[[i]]
    all_cor_pot <- pot_corr_list[[i]]
    num_gene_exp <- length(all_enh_pot)
    index_list <- as.list(1:length(all_enh_pot))
    use_pos<- use_pos_list[[i]]
    #07252019
    enh_act <- enh_act_list[[i]]
    gene_act_obs <- gene_act[[i]]
    
    all_group_tf <- list()
    #(11/15) all_group_dist <- list()
    all_gene_tf <- list()
    all_group_corr <- list()
    all_group_non_corr <- list()
    all_group_D_scalar <- list()
    for (gene_num in 1:num_gene_exp) {
      gene_name <- all_gene_name[gene_num]
      pos_in_whole <- which(gene_infor_matrix[,1] == gene_name)
      which_group <- group_infor[pos_in_whole]
      group_tf <- group_tf_mean_data[[which_group]] ##(11/17)
      #(11/15) group_dist <- group_tf_dist[[which_group]]
      gene_tf <- TF_GE[gene_num, ]
      group_corr <- group_corr_dist[[which_group]]
      group_non_corr <- group_non_corr_dist[[which_group]]
      all_group_D_scalar[[gene_num]] <- calculate_D(entropy(group_tf+mean_group_tf))
      all_group_tf[[gene_num]] <- group_tf
      #(11/15) all_group_dist[[gene_num]] <- group_dist
      all_gene_tf[[gene_num]] <- gene_tf
      all_group_corr[[gene_num]] <- group_corr
      all_group_non_corr[[gene_num]] <- group_non_corr 
    }
    
    ## major function(apply to each gene)
    update_gene <- function(pot_enh, inter_enh, distance ,corr, group_tf, group_corr, group_non_corr, one_TF_GE, enh_tf,D_scalar,enh_act,gene_act_obs,ind,all_tf_count) {
      #print(ind)
      ## function: to calculate the log-probability a enhancer blongs to a gene(fixed)
      prob_link<- function(enh_num) {
        
        #one_TF_GE <- TF_GE[gene_num, ]
        
        # pot_enh <- all_enh_pot[[gene_num]]
        # inter_enh <- all_enh_inter[[gene_num]]
        # corr <- all_cor_pot[[gene_num]]
        
        pos_this_enh <- which(pot_enh == enh_num)
        corr_obs <- corr[pos_this_enh]
        dist_obs <- distance[pos_this_enh]
        #0902
        #TF_ENH <- enh_tf[which(use_pos==enh_num), ]
        
        #07252019: if activity is 0, return 0
        enh_act_obs <- enh_act[pos_this_enh]
        act <- enh_act_obs*gene_act_obs
        tmp_tf_count <- all_tf_count[pos_this_enh]
        
        #changed on 0905, to set a prob for no motif one
        if(act==0){
          
          return(rep(-1e6,3))
          
        }else{
          
          TF_ENH <- enh_tf[which(use_pos==enh_num), ]
          
        }
        
        if (enh_num %in% inter_enh) {
          new_obs <- one_TF_GE
          new_obs <- new_obs/tf_sum #(11/21)
          new_obs <- replace_inf(new_obs) #(11/21)
          new_obs <- replace_na(new_obs) #(11/21)
          #prob_one <- get_log_prob(group_dist, new_obs)
          #prob_one <- -min(KL.Dirichlet(new_obs, group_tf,1,1), KL.Dirichlet(group_tf, new_obs,1,1))
          #prob_one <- -KLD(new_obs,group_tf)[[5]]
          
          #0905 
          if(tmp_tf_count!=0){
            prob_one <- -KL(new_obs,group_tf)
            prob_one <- prob_one * D * D_scalar * d_m}else{
              
              prob_one <- log(0.46)
              
            }
          
          
          #(01/24) corr_infor <- as.matrix(corr_inter[-which(corr_inter == corr_obs)])
          #(01/24) prob_two <- C * log_prob_occur(corr_infor, corr_obs)
          prob_two <- C * get_log_prob(group_corr, corr_obs)
          prob_three <- B * get_log_prob(inter_distance, dist_obs)
        } else {
          new_obs <- add_one_mean(one_TF_GE, num_int_enh, TF_ENH)
          new_obs <- new_obs/tf_sum #(11/21)
          new_obs <- replace_inf(new_obs) #(11/21)
          new_obs <- replace_na(new_obs) #(11/21)
          #prob_one <- get_log_prob(group_dist, new_obs)
          #prob_one <- -min(KL.Dirichlet(new_obs, group_tf,1,1), KL.Dirichlet(group_tf, new_obs,1,1))
          #prob_one <- -KLD(new_obs,group_tf)[[5]]
          #0905 
          if(tmp_tf_count!=0){
            prob_one <- -KL(new_obs,group_tf)
            prob_one <- prob_one * D * D_scalar * d_m}else{
              
              prob_one <- log(0.53)
              
            }
          #(01/24) corr_infor <- as.matrix(corr_inter)
          #(01/24) prob_two <- C * log_prob_occur(corr_infor, corr_obs)
          prob_two <- C * get_log_prob(group_corr, corr_obs)
          prob_three <- B * get_log_prob(inter_distance, dist_obs)
          #prior <- log(length(unlist(inter_enh_list)))
          
        }
        
        return(c(prob_one, prob_two, prob_three))
        
      }
      
      ## function: to calculate the log-probability a enhancer doesn't blongs to a gene(fixed)
      prob_no_link<- function(enh_num) {
        
        #one_TF_GE <- TF_GE[gene_num, ]
        
        # pot_enh <- all_enh_pot[[gene_num]]
        # inter_enh <- all_enh_inter[[gene_num]]
        # corr <- all_cor_pot[[gene_num]]
        
        pos_this_enh <- which(pot_enh == enh_num)
        corr_obs <- corr[pos_this_enh]
        dist_obs <- distance[pos_this_enh]
        
        #0902
        #TF_ENH <- enh_tf[use_pos[which(use_pos==enh_num)], ]
        
        #07252019: if activity is 0, return 0
        enh_act_obs <- enh_act[pos_this_enh]
        act <- enh_act_obs*gene_act_obs
        tmp_tf_count <- all_tf_count[pos_this_enh]
        if(act==0){
          
          return(rep(0,3))
          
        }else{
          
          TF_ENH <- enh_tf[which(use_pos==enh_num), ]
          
        }
        
        if (enh_num %in% inter_enh) {
          new_obs <- remove_one_mean(one_TF_GE, num_int_enh, TF_ENH)
          new_obs <- new_obs/tf_sum #(11/21)
          new_obs <- replace_inf(new_obs) #(11/21)
          new_obs <- replace_na(new_obs) #(11/21)
          #prob_one <- get_log_prob(group_dist, new_obs)
          #prob_one <- -min(KL.Dirichlet(new_obs, group_tf,1,1), KL.Dirichlet(group_tf, new_obs,1,1))
          #prob_one <- -KLD(new_obs,group_tf)[[5]]
          
          #0905 
          if(tmp_tf_count!=0){
            prob_one <- -KL(new_obs,group_tf)
            prob_one <- prob_one * D * D_scalar * d_m}else{
              
              prob_one <- log(0.57)
              
            }
          
          #(01/24) corr_infor <- as.matrix(corr_no_inter)
          #(01/24) prob_two <- C * log_prob_occur(corr_infor, corr_obs)
          prob_two <- C * get_log_prob(group_non_corr, corr_obs)
          prob_three <- B * get_log_prob(back_distance, dist_obs)
        } else {
          new_obs <- one_TF_GE 
          new_obs <- new_obs/tf_sum #(11/21)
          new_obs <- replace_inf(new_obs) #(11/21)
          new_obs <- replace_na(new_obs) #(11/21)
          #prob_one <- get_log_prob(group_dist, new_obs)
          #prob_one <- -min(KL.Dirichlet(new_obs, group_tf,1,1), KL.Dirichlet(group_tf, new_obs,1,1))
          #prob_one <- -KLD(new_obs,group_tf)[[5]]
          
          #0905 
          if(tmp_tf_count!=0){
            prob_one <- -KL(new_obs,group_tf)
            prob_one <- prob_one * D * D_scalar * d_m}else{
              
              prob_one <- log(0.54)
              
            }
          
          #(01/24) corr_infor <- as.matrix(corr_no_inter[-which(corr_no_inter == corr_obs)])
          #(01/24) prob_two <- C * log_prob_occur(corr_infor, corr_obs)
          prob_two <- C * get_log_prob(group_non_corr, corr_obs)
          prob_three <- B * get_log_prob(back_distance, dist_obs)
          #prior <- log(length(unlist(non_inter_corr_list)))
        }
        
        return(c(prob_one, prob_two, prob_three))
        
      }
      
      ## Function: To update enhancer-gene interaction
      update_enh <- function(enh_num) {
        # vector_link<- prob_link(enh_num)
        # vector_no_link<- prob_no_link(enh_num)
        # log_probs_matrix<- rbind(vector_link, vector_no_link)
        
        vector_link<- sum(prob_link(enh_num))
        vector_no_link<- sum(prob_no_link(enh_num))
        log_probs_matrix<- c(vector_link, vector_no_link)
        
        
        # probs_sep<- normalize_log_prob(log_probs_matrix)
        # probs<- as.vector(unlist(apply(probs_sep, 1, FUN = prod)))
        probs <- as.vector(normalize_log_prob(log_probs_matrix))
        
        probs <- normalize_to_one(probs)
        
        new_result<- sample(c(1,0), size = 1, prob = c(probs[1], 1-probs[1]))
        
        #03/27/2019 need posterior probability
        #if (new_result == 1) {
        #  return(enh_num)
        #} else {
        #  return(0)
        #}
        
        return(probs)
        
      } 
      update_enh <- Vectorize(update_enh)
      
      num_int_enh<- length(inter_enh)
      
      pos_inter<- which(pot_enh %in% inter_enh)
      corr_inter<- corr[pos_inter]
      corr_no_inter<- corr[-pos_inter]
      
      new_inter_enh <- update_enh(pot_enh)
      
      #03272019 sample label from posterior probability
      
      post_p <- new_inter_enh
      
      new_inter_enh <- apply(post_p,2,function(x){
        
        sample(c(1,0),size=1,prob=c(x[1],x[2]))
        
      })
      
      
      #new_inter_enh <- remove_zero(new_inter_enh)
      
      #03272019
      new_inter_enh <- pot_enh[which(new_inter_enh>0)]
      if(length(new_inter_enh)==0){new_inter_enh <- NULL}
      
      
      #return(as.vector(new_inter_enh))
      
      #03272019 return a list, first entry is label, second entry is post_p
      
      return(list(new_inter_enh,post_p))
      
    }
    
    
    new_all_enh_inter <- mcmapply(FUN = update_gene, pot_enh = all_enh_pot, inter_enh = all_enh_inter, distance = all_dist_pot, corr = all_cor_pot, 
                                  group_tf = all_group_tf, group_corr = all_group_corr, group_non_corr = all_group_non_corr, 
                                  one_TF_GE = all_gene_tf,D_scalar = all_group_D_scalar,gene_act_obs = gene_act_obs,enh_act = enh_act, ind=index_list,all_tf_count=all_tf_count,MoreArgs=list(enh_tf), mc.cores = num_cores)
    
    ## (01/24)
    result_type <- unlist(lapply(new_all_enh_inter, class)) == "character"
    ## (01/18)
    r <- 1
    while (sum(result_type)>0 & r<10000) {
      print(paste("we got an error", r))
      new_all_enh_inter <- mcmapply(FUN = update_gene, pot_enh = all_enh_pot, inter_enh = all_enh_inter, corr = all_cor_pot, 
                                    group_tf = all_group_tf, group_corr = all_group_corr, group_non_corr = all_group_non_corr, 
                                    one_TF_GE = all_gene_tf,D_scalar = all_group_D_scalar,gene_act = gene_act,enh_act = enh_act_list,MoreArgs=list(enh_tf), mc.cores = num_cores, mc.preschedule=F)
      result_type <- unlist(lapply(new_all_enh_inter, class)) == "character"
      r <- r +1
    }
    # end_time <- Sys.time()
    
    #03272019 
    post_p[[i]] <- lapply(seq(2,length(new_all_enh_inter),2),function(x){new_all_enh_inter[[x]]})
    new_all_enh_inter <- lapply(seq(1,length(new_all_enh_inter),2),function(x){new_all_enh_inter[[x]]})
    
    
    new_bin_link <- mapply(FUN=in_or_not, x=all_enh_pot, y=new_all_enh_inter)
    ############################################### update interaction #######################################
    inter_enh_list[[i]] <- new_all_enh_inter
    bin_link_list[[i]] <- new_bin_link
    
    
    
  } ## chr-loop ends 
  print(paste("all enhancer-gene interaction update done"))
  enh_end <- Sys.time() ## update enhancer-gene interaction end time
  # print(paste("all enhancer-gene interaction update done in", enh_end-enh_start))
  
  # save(
  #   inter_enh_list,
  #   post_p,
  #   file = paste(dir_path, "/inter_enh_list.Rdata",sep="") 
  # )
  
  
  ######################## To update gene group membership ###########################
  
  ## update TF_GE_list for all genes in all chromosomes (no normalization)
  all_gene_group_mean_tf <- all_gene_group_mean_tf <- apply(all_TF_GE,2,function(x){mean(x,na.rm=T)})
  new_TF_GE_list <- list()
  for (i in 1:length) {
    all_enh_inter_exp <- inter_enh_list[[i]]
    all_enh_pot_exp <- pot_enh_list[[i]]
    enh_tf <- enh_tf_list[[i]]
    gene_name_all_exp <- gene_name_list[[i]]
    
    num_gene_exp <- length(gene_name_all_exp)
    num_of_tf <- dim(enh_tf)[2]
    
    TF_GE <- matrix(0,num_gene_exp,num_of_tf)
    for (g in 1:num_gene_exp) {
      inter_enh_pos <- match(all_enh_inter_exp[[g]],use_pos_list[[i]])
      
      #0905 for enhancer with no TF, will get a NA
      sub_enh_pos <- inter_enh_pos[which(!is.na(inter_enh_pos))]
      
      if (length(sub_enh_pos) != 0) {
        sub_enh_tf <- matrix(enh_tf[sub_enh_pos, ], length(sub_enh_pos), num_of_tf)
        TF_GE[g, ] <- as.vector(unlist(apply(sub_enh_tf, 2, mean)))
      } else {
        pot_enh_pos <- match(all_enh_pot_exp[[g]],use_pos_list[[i]])
        pot_enh_pos <- pot_enh_pos[which(!is.na(pot_enh_pos))]
        sub_enh_tf <- matrix(enh_tf[pot_enh_pos, ], length(pot_enh_pos), num_of_tf)
        TF_GE[g, ] <- as.vector(unlist(apply(sub_enh_tf, 2, mean)))
        if(length(sub_enh_tf)==0){TF_GE[g,] <- all_gene_group_mean_tf}
      }
      
    }
    
    rownames(TF_GE) <- gene_name_all_exp
    colnames(TF_GE) <- colnames(enh_tf)
    new_TF_GE_list[[i]] <- TF_GE
  }
  
  
  
  ############################################ remember previous result #################################
  pre_TF_GE_list <- TF_GE_list
  
  pre_all_TF_GE <- all_TF_GE
  
  pre_tf_sum <- tf_sum
  
  pre_count_TF_GE <- count_TF_GE
  
  pre_all_gene_tf <- all_gene_tf
  
  pre_inter_corr_list <- inter_corr_list
  
  pre_non_inter_corr_list <- non_inter_corr_list
  
  # pre_all_gene_corr  <- all_gene_corr 
  ########################################### update previous result #####################################
  TF_GE_list <- new_TF_GE_list
  
  
  ## prepare all_TF_GE matrix and list(noramlized by column/TF)
  new_all_TF_GE <- do.call(rbind,TF_GE_list)
  all_TF_GE <- new_all_TF_GE
  tf_sum <- apply(all_TF_GE, 2, sum) #(11/21)
  count_TF_GE <- all_TF_GE #(11/21)
  all_TF_GE <- apply(all_TF_GE, 2, prop)
  all_TF_GE <- apply(all_TF_GE, 2, replace_na)
  all_gene_tf <- trans_to_list(all_TF_GE,1)
  
  
  ## update inter_corr_list and non_inter_corr_list (seperated by chromosome)
  for (i in 1:length) {
    all_enh_inter_exp <- inter_enh_list[[i]]
    all_enh_pot_exp <- pot_enh_list[[i]]
    all_cor_pot <- pot_corr_list[[i]]
    num_gene_exp <- length(all_enh_pot_exp)
    all_act <- pot_act_list[[i]]
    all_gene_act <- gene_act[[i]]
    all_tf_count <- tf_count_list[[i]]
    all_inter_corr <- list()
    all_non_inter_corr <- list()
    for (g in 1:num_gene_exp) {
      inter <- all_enh_inter_exp[[g]]
      pot <- all_enh_pot_exp[[g]]
      corr <- all_cor_pot[[g]]
      act_temp <- all_act[[g]]
      tf_count_ind <- all_tf_count[[g]]
      tmp_gene_act <- all_gene_act[[g]]
      if (length(inter) > 0 ) {
        pos <- which(pot %in% inter)
        act <- which(act_temp*tmp_gene_act>0)
        tf_ind <- which(tf_count_ind>0)
        
        #0905 relax to active enhancers
        non_link_ind <- setdiff(act,pos)
        all_inter_corr[[g]] <- corr[pos]
        all_non_inter_corr[[g]] <- c(corr[non_link_ind],logical(0))
      } else {
        #all_inter_corr[[g]] <- rep(0, pcount_corr)
        act <- which(act_temp*tmp_gene_act>0)
        tf_ind <- which(tf_count_ind>0)
        all_inter_corr[[g]] <- rep(mean(corr), pcount_corr)
        all_non_inter_corr[[g]] <- c(corr[act],logical(0))
      }
    }
    inter_corr_list[[i]] <- all_inter_corr
    non_inter_corr_list[[i]] <- all_non_inter_corr
  }
  
  
  ## prepare all_inter_corr list (whole geno)
  all_gene_corr <- list()
  for (gene_num in 1:dim(gene_infor_matrix)[1]) {
    which_chr <- gene_infor_matrix[gene_num,2]
    pos_in_chr <- gene_infor_matrix[gene_num,3]
    corr_obs <- inter_corr_list[[which_chr]][[pos_in_chr]]
    all_gene_corr[[gene_num]] <- corr_obs
  }
  
  
  
  
  
  ## To update gene group membership
  
  # start_time <- Sys.time()
  
  # ## use naive hierarchical cluster to update group membership
  # data <- all_TF_GE
  # data <- data[,-as.vector(which(apply(data,2,sum) == 0))]
  # phm <- pheatmap(data,color = colorRampPalette(c("snow", "red"))(50),clustering_distance_row = "correlation", clustering_distance_col = "correlation", fontsize=9, fontsize_row=1, fontsize_col=3)
  # tree_gene <- phm$tree_row
  # tree_tf <- phm$tree_col
  # new_group_infor <- as.vector(cutree(tree_gene, k = num_of_group))
  # 
  
  ## function: calculate log-probability a gene belongs to a group
  ## return the probabilities one gene belongs to all groups respectively
  update_group <- function(gene_tf, gene_corr, gene_group, gene_name) {
    
    #this_group_tf <- group_tf_data[[gene_group]]
    #pos_this_gene_in_group <- which(rownames(this_group_tf) == gene_name)
    #this_group_tf <- this_group_tf[-pos_this_gene_in_group,]
    #new_one_tf_dist <- get_log_dist(this_group_tf)
    
    #this_group_corr <- group_corr_data[[gene_group]]
    #pos_this_gene_corr <- NULL
    #for (l in 1:length(gene_corr)) {
    #  pos <- which(this_group_corr == gene_corr[l])[1]
    #  pos_this_gene_corr <- c(pos_this_gene_corr, pos)
    #}
    #this_group_corr <- this_group_corr[-pos_this_gene_corr]
    #new_one_corr_dist <- get_log_dist(this_group_corr)
    
    ##new_tf_dist <- group_tf_dist
    #new_tf_dist[[gene_group]] <- new_one_tf_dist 
    ##new_corr_dist <- group_corr_dist
    #new_corr_dist[[gene_group]] <- new_one_corr_dist
    
    
    tf_mean_data <- group_tf_mean_data
    prob_belong <- function(one_group_tf) {
      
      #prob_left <- get_log_prob(one_group_tf, gene_tf)
      #(11/15)prob_left <- -min(KL.Dirichlet(gene_tf, one_group_tf,1,1), KL.Dirichlet(one_group_tf, gene_tf,1,1))
      #prob_left <- -KLD(gene_tf, one_group_tf)[[5]]
      prob_left <- -KL(one_group_tf, gene_tf)
      prob_left <- prob_left * DD
      
      # (02/20)
      # prob_corr <- NULL
      # for (l in 1:length(gene_corr)) {
      #   one_corr <- gene_corr[l]
      #   one_prob <- get_log_prob(one_group_corr, one_corr)
      #   prob_corr <- c(prob_corr, one_prob)
      # }
      # prob_right <- sum(prob_corr)
      
      #(10/09) return(c(prob_left, prob_right))
      
      # (02/20) return(sum(prob_left, prob_right))
      return(prob_left)
    }
    
    # (02/20) one_gene_log_probs <- mapply(FUN = prob_belong, one_group_tf = tf_mean_data, one_group_corr = group_corr_dist)
    one_gene_log_probs <- mapply(FUN = prob_belong, one_group_tf = tf_mean_data)
    
    
    #(10/09) one_gene_log_probs <- t(one_gene_log_probs)
    one_gene_log_probs <- as.matrix(one_gene_log_probs, num_of_group,1)
    
    
    #(10/09) probs_sep <- normalize_log_prob(one_gene_log_probs)
    # probs <- as.vector(unlist(apply(probs_sep, 1, FUN = prod)))
    probs <- as.vector(normalize_log_prob(one_gene_log_probs))
    
    
    one_gene_probs <- normalize_to_one(probs)
    #one_gene_probs <- probs
    
    return(one_gene_probs)
    
  }
  
  all_groups_probs <- mcmapply(FUN = update_group, gene_tf = all_gene_tf, gene_corr = all_gene_corr, gene_group = group_infor_list, gene_name = all_gene_name_list, mc.cores = num_cores)
  
  new_group_infor <- apply(all_groups_probs, 2, samp_one)
  
  #all_groups_probs_list <- trans_to_list(all_groups_probs,2)
  #new_group_infor <- unlist(mclapply(all_groups_probs_list, samp_one, mc.cores = num_cores))
  
  #################################################### remember old gene group data #######################  
  
  pre_group_infor <- group_infor
  
  pre_group_tf_mean_data <- group_tf_mean_data
  
  pre_group_corr_data <- group_corr_data
  
  pre_group_non_corr_data <- group_non_corr_data
  
  pre_group_corr_dist <- group_corr_dist
  
  pre_group_non_corr_dist <- group_non_corr_dist
  
  pre_group_tf_data <- group_tf_data
  
  ################################################### update gene group data ##############################
  
  
  group_infor <- new_group_infor
  
  ## prepare gene group list
  group_infor_list <- list()
  for (g in 1:length(group_infor) ) {
    group_infor_list[[g]] <- group_infor[g]
  }
  
  
  ## prepare TF data  for each gene group
  mean_gene_tf <- apply(all_TF_GE, 2, mean)
  group_tf_data <- list()
  for (k in 1:num_of_group) {
    
    gene_pos <- which(group_infor == k)
    if (length(gene_pos) > 0) {
      
      if (length(gene_pos) > 1) {
        tf_for_this_group <- all_TF_GE[gene_pos, ]
      } else {
        tf_for_this_group <- rbind(all_TF_GE[gene_pos, ], mean_gene_tf)
      }
      group_tf_data[[k]] <- tf_for_this_group
      
    } else {
      #group_tf_data[[k]] <- matrix(0, pcount_gene, dim(all_TF_GE)[2])
      group_tf_data[[k]] <- matrix(rep(mean_gene_tf, pcount_gene), pcount_gene, dim(all_TF_GE)[2], byrow = T)
    }
  }
  
  
  
  group_tf_mean_data <- lapply(group_tf_data, take_mean_in_list)
  
  
  # (11/15)
  # ## prepare TF distribution for each gene group
  # group_tf_dist <- list()
  # for (k in 1:num_of_group) {
  #   group_tf<- group_tf_data[[k]]
  #   group_tf_dist[[k]] <- get_log_dist(group_tf)
  # }
  
  
  ## prepare correlation data for each gene group
  mean_gene_corr <- mean(unlist(inter_corr_list))
  group_corr_data <- list()
  for (k in 1:num_of_group) {
    gene_pos <- which(group_infor == k)
    
    if (length(gene_pos) != 0 ) {
      
      corr_one_group <- NULL
      for (l in gene_pos) {
        which_chr <- gene_infor_matrix[l,2]
        pos_in_chr <- gene_infor_matrix[l,3]
        corr_one_gene_inter <- inter_corr_list[[which_chr]][[pos_in_chr]]
        corr_one_group <- c(corr_one_group, corr_one_gene_inter)
      }
      if (length(corr_one_group) == 1) {
        corr_one_group <- c(corr_one_group, mean_gene_corr)
      }
      group_corr_data[[k]] <- as.matrix(corr_one_group)
      
    } else {
      #group_corr_data[[k]] <- rep(0, pcount_corr)
      group_corr_data[[k]] <- as.matrix(rep(mean_gene_corr, pcount_corr))
    }
    
  }
  
  ## prepare non-intercted correlation data for each gene group
  mean_non_gene_corr <- mean(unlist(non_inter_corr_list))
  group_non_corr_data <- list()
  for (k in 1:num_of_group) {
    gene_pos <- which(group_infor == k)
    
    if (length(gene_pos) != 0 ) {
      
      corr_one_group <- NULL
      for (l in gene_pos) {
        which_chr <- gene_infor_matrix[l,2]
        pos_in_chr <- gene_infor_matrix[l,3]
        corr_one_gene_non_inter <- non_inter_corr_list[[which_chr]][[pos_in_chr]]
        #if(any(is.na(corr_one_gene_non_inter))){break}
        corr_one_group <- c(corr_one_group, corr_one_gene_non_inter)
      }
      if (length(corr_one_group) == 1) {
        corr_one_group <- c(corr_one_group, mean_non_gene_corr)
      }
      group_non_corr_data[[k]] <- as.matrix(corr_one_group)
      
    } else {
      #group_corr_data[[k]] <- rep(0, pcount_corr)
      group_non_corr_data[[k]] <- as.matrix(rep(mean_non_gene_corr, pcount_corr))
    }
    
  }
  
  
  
  ## prepare correlation distribution for each gene group
  group_corr_dist <- list()
  data <- unlist(inter_corr_list)
  d <- get_log_dist(data)
  for (k in 1:num_of_group) {
    one_group_dist <- d
    group_corr_dist[[k]] <- one_group_dist
  }
  
  ## prepare non-interacted correlation distribution for each gene group
  group_non_corr_dist <- list()
  data <- unlist(non_inter_corr_list)
  d <- get_log_dist(data)
  for (k in 1:num_of_group) {
    one_group_dist <- d
    group_non_corr_dist[[k]] <- one_group_dist
  }
  
  
  ################################################# simulated annealling ##################################
  
  #calculate engergy 
  e <- calculate_energy(pre_group_tf_mean_data,pre_group_corr_data,pre_group_non_corr_data,K)
  
  if(n > 1){
    e_star <- calculate_energy(group_tf_mean_data,group_corr_data,group_non_corr_data,K)
  }else{
    
    e_star <- e
    
  }
  
  #calculate accept probability
  
  p <- accept_prob(e,e_star,Temperature=Temperature)
  
  #decide reject or accept
  acc <- decision(p)
  print(acc)
  print(e_star)
  print(e)
  if(acc == 1){
    #if accept nothing need to update, recalculate T
    Temperature <- Temperature*(100-n)/100
    if(n==1){
      D_pre <- D
    }
    D <- D
    print(paste("accept, new D = ",D,"former D = ",D_pre,"accept_probability:",p))
    
  }else{
    #if reject, go to previous result
    inter_enh_list <- pre_enh_list
    
    inter_corr_list <- pre_inter_corr_list
    
    post_p <- pre_post_p
    
    bin_link_list <- pre_bin_link_list
    
    TF_GE_list <- pre_TF_GE_list
    
    all_TF_GE <- pre_all_TF_GE
    
    tf_sum <- pre_tf_sum
    
    count_TF_GE <- pre_count_TF_GE
    
    group_infor <- pre_group_infor
    
    group_tf_mean_data <- pre_group_tf_mean_data
    
    group_tf_data <- pre_group_tf_data
    
    group_corr_data <- pre_group_corr_data
    
    group_non_corr_data <- pre_group_non_corr_data
    
    group_corr_dist <- pre_group_corr_dist
    
    group_non_corr_dist <- pre_group_non_corr_dist
    
    all_gene_tf <- pre_all_gene_tf
    
    non_inter_corr_list <- pre_non_inter_corr_list
    
    # all_gene_corr  <- pre_all_gene_corr
    
    #0521 if reject, temperature don't change
    #Temperature <- Temperature*(100-n)/100
    Temperature <- Temperature
    print(paste("reject, new D = ",D_pre,"sampled D = ",D,"accept_probability:",p))
    
    D <- D_pre
  }
  
  #sample a new D
  D_pre <- D
  
  D <- rtruncnorm(1,D_range[1],D_range[2],D_pre,Temperature)
  
  if(abs(D-D_pre)<1e-4) break
  ######################## save data for this iteration ###########################
  iter_end_time <- Sys.time()
  print(paste("iteration", n, "finished in: "), sep="")
  print(iter_end_time - iter_start_time)
  
  save(
    D_pre,
    DD,
    C,
    acc,
    
    num_of_group,
    gene_infor_matrix,
    
    inter_corr_list,
    inter_enh_list,
    non_inter_corr_list,
    post_p,
    bin_link_list,
    pre_bin_link_list,
    
    TF_GE_list,
    all_TF_GE,
    tf_sum,
    count_TF_GE,
    
    all_groups_probs,
    group_infor,
    group_tf_data,
    pre_group_tf_data,
    group_corr_data,
    group_non_corr_data,
    
    file = file_name_iter[n]
    
  )
  
}## iteration loop ends
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##
## iteration ends here
##
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




## update some group-related data after all iterations
## prepare gene group list
group_infor_list <- list()
for (g in 1:length(group_infor) ) {
  group_infor_list[[g]] <- group_infor[g]
}


## prepare TF data  for each gene group
mean_gene_tf <- apply(all_TF_GE, 2, mean)
group_tf_data <- list()
for (k in 1:num_of_group) {
  
  gene_pos <- which(group_infor == k)
  if (length(gene_pos) > 0) {
    
    if (length(gene_pos) > 1) {
      tf_for_this_group <- all_TF_GE[gene_pos, ]
    } else {
      tf_for_this_group <- rbind(all_TF_GE[gene_pos, ], mean_gene_tf)
    }
    group_tf_data[[k]] <- tf_for_this_group
    
  } else {
    #group_tf_data[[k]] <- matrix(0, pcount_gene, dim(all_TF_GE)[2])
    group_tf_data[[k]] <- matrix(rep(mean_gene_tf, pcount_gene), pcount_gene, dim(all_TF_GE)[2], byrow = T)
  }
}


group_tf_mean_data <- lapply(group_tf_data, take_mean_in_list)


## prepare correlation data for each gene group
mean_gene_corr <- mean(unlist(inter_corr_list))
group_corr_data <- list()
for (k in 1:num_of_group) {
  gene_pos <- which(group_infor == k)
  
  if (length(gene_pos) != 0 ) {
    
    corr_one_group <- NULL
    for (l in gene_pos) {
      which_chr <- gene_infor_matrix[l,2]
      pos_in_chr <- gene_infor_matrix[l,3]
      corr_one_gene_inter <- inter_corr_list[[which_chr]][[pos_in_chr]]
      corr_one_group <- c(corr_one_group, corr_one_gene_inter)
    }
    
    if (length(corr_one_group) == 1) {
      corr_one_group <- c(corr_one_group, mean_gene_corr)
    }
    
    group_corr_data[[k]] <- as.matrix(corr_one_group)
    
  } else {
    #group_corr_data[[k]] <- rep(0, pcount_corr)
    group_corr_data[[k]] <- as.matrix(rep(mean_gene_corr, pcount_corr))
  }
  
}


save(
  D,
  DD,
  C,
  
  num_of_group,
  gene_infor_matrix,
  
  inter_corr_list,
  inter_enh_list,
  non_inter_corr_list,
  bin_link_list,
  pre_bin_link_list,
  
  TF_GE_list,
  all_TF_GE,
  tf_sum,
  count_TF_GE,
  
  all_groups_probs,
  group_infor,
  pre_group_tf_data,
  group_tf_data,
  group_corr_data,
  group_non_corr_data,
  
  file = file_name_iter[num_iter]
  
)


all_end_time <- Sys.time()
print(paste("all iterations finished in: "))
print(all_end_time - all_start_time)

