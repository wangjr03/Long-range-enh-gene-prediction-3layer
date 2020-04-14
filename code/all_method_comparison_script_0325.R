args <- commandArgs(T)

ripple_link <- args[1]
  
focs_link <- args[2]

ernst_link <- args[3]
  
targetfinder_link <- args[4]
  
  
jeme_link <- args[5]

im_pet_link <- args[6]
  
info_frame_link <- args[7]

link_frame_link <- args[8]

gd <- args[9]

figure_name <- args[10]

p_w <- as.numeric(args[11])

###
load(info_frame_link)

#info_frame <- chunk.merge


link_frame <- read.table(link_frame_link)

#link_frame <- read.table("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/Link_frame/link_frame.KejiZhao_hg19_CD4+T.bed.bed")


#gd <- 6


#

tmp.pth <- getwd()

setwd(gd)


n <- length(dir())

ind <- list()
for(i in 1:n){
  
  ind[[i]] <- data.table::fread(paste0(i,".txt"))
  
}
ind <- do.call(rbind,ind)
ind <- as.data.frame(ind)

info_frame$new_col <- ind[,1]

gd <- dim(info_frame)[2]

setwd(tmp.pth)




promoter_frame <- read.table("../promoter_activity/promoter_frame.txt")
promoter_frame <- promoter_frame[,c(2,3,4,1)]
#compare ripple
print("comparing Ripple")
if(ripple_link!="NA"){
  
  #ripple <- read.table("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/RIPPLE/predicitons/Gm12878.txt")
  ripple <- read.table(ripple_link)
  
  
  ripple_enh <- as.data.frame(do.call(rbind,strsplit(as.character(ripple[,1]),"_")))
  
  ripple_enh <- data.frame(ripple_enh)
  
  ripple_enh[,2] <- as.numeric(as.character(ripple_enh[,2]))
  
  ripple_enh[,3] <- as.numeric(as.character(ripple_enh[,3]))
  gene_annotation <- read.table("../gene_annotation/gene_annotation_V19.txt")
  ripple_gene_name <- c()
  
  for(i in 1:dim(ripple)[1]){
    
    name <-  as.character(gene_annotation[which(gene_annotation[,6]==as.character(ripple[i,3])),7])
    
    if(length(name)==0){
      
      name <- "NA"
      
    }else{
    
    name <- name[1]
    
    }
    
    ripple_gene_name <- c(ripple_gene_name,name)
    
  }
  ripple <- cbind(ripple,ripple_enh)
  
  
  ripple$gene_name <- ripple_gene_name
  
  #alt_id <- grep(as.character(ripple$'V3'[which((ripple$gene_name=="NA"))]),gene_annotation[,7])
  
  #ripple$gene_name[which((ripple$gene_name=="NA"))] <- as.character(ripple$'V3'[which((ripple$gene_name=="NA"))])
  
  ensg_name_loc <- grep("ENSG",ripple$'V3')
  
  ensg_name <- as.character(ripple$'V3'[ensg_name_loc])
  
  full_name <- sapply(ensg_name,function(x){
    
    tmp_id <- grep(x,gene_annotation[,7])
    
    if(length(tmp_id)==0){"NA"}else{
    
    as.character(gene_annotation[tmp_id[1],7])
    
    }
    
  })
  
  ripple$gene_name[ensg_name_loc] <- full_name
  
  names(ripple) <- c("enh_coord","ripple_score","gene","chr","start","stop","ensg_gene")
  
  # 
  # ripple <- cbind(ripple,ripple_enh)
  # 
  # ripple$V7 <- ripple_gene_name
  # 
  
  ripple[,5] <- ripple[,5]+150
  
  ripple[,6] <- ripple[,6]-150
  
  
  
  match_ripple <- function(x){
    
    gene_loc <- gene_annotation[which(gene_annotation[,7]==ripple[x,7]),]
    
    if(dim(gene_loc)[1]==0){return(FALSE)}
    
    if(gene_loc[4]==1){
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,2]-1000,stop=gene_loc[,2]+1000)
      
    }else{
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
      
    }
    
    enh_loc <- ripple[x,c(4,5,6)]
    
    ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
                link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))
    any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
          link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
    
    return(ind)
    
  }
  
  match_ripple <- Vectorize(match_ripple)
  
  ripple_ind <- match_ripple(1:dim(ripple)[1])
  tss <- list()
  #not control distance, ripple Thymus GM12878 0.015
  
  for(x in 1:dim(ripple)[1]){
    
    
    gene_loc <- gene_annotation[which(gene_annotation[,7]==ripple[x,7]),]
    
    if(dim(gene_loc)[1]==0){tss[[x]] <- NA;next}
    
    gene_loc <- gene_loc[1,]
    
    if(gene_loc[4]==1){
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,2]-1000,stop=gene_loc[,2]+1000)
      
    }else{
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
      
    }
    
    tss[[x]] <- gene_loc
    
  }
  
  tss <- do.call(rbind,tss)
  
  ripple <- cbind(ripple,tss)
  
  my_pair <- list()
  
  
  ripple_d <- abs(ripple[,9]-ripple[,5])
  
  my_pair <- list()
  
  #my_pair.2 <- list()
  
  sub_info_frame <- subset(info_frame,info_frame$dist<1e6&info_frame$act==1)
  
  
  for(i in 1:100){
    #print(i)
    d_s <- (i-1)*10000+1
    
    d_e <- i*10000
    
    count <- sum(ripple_d<d_e&ripple_d>d_s,na.rm = T)
    
    my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
    
    n <- dim(my_pool)[1]
    #my_pred <- my_pool[order(my_pool$adj_p,decreasing = F),][1:count,]
    
    my_pred <- my_pool[order(my_pool$corrected_p,decreasing = T),][1:min(n,count),]
    
    my_pair[[i]] <- my_pred
    #my_pair.2[[i]] <- my_pred.2
    
  }
  
  my_pair <- data.table::rbindlist(my_pair)
  
  #my_pair.2 <- data.table::rbindlist(my_pair.2)
  
  my_pair <- as.data.frame(my_pair)
  
  #my_pair.2 <- as.data.frame(my_pair.2)
  
  ripple_my_frac <- mean(my_pair[,gd])
  
  #ripple_my_frac.2 <- mean(my_pair.2[,gd])
}else{
  
  ripple_ind <- 0
  
  ripple_my_frac <- 1
  
  
}

#compare FOCS
print("comparing FOCS")
#focs <- read.table("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/FOCS/FOCS-add-license-1/predictions/encode_interactions.txt",sep="\t",header=T)
if(focs_link!="NA"){

focs <- read.table(focs_link,sep="\t",header=T)

#resize focs enhancer and promoter size

p_center <- floor((focs$p_start + focs$p_end)/2)

e_center <- floor((focs$e_start + focs$e_end)/2)

focs[,2] <- p_center-1000

focs[,3] <- p_center+1000

focs[,6] <- e_center-500

focs[,7] <- e_center+500

focs[,1] <- as.character(focs[,1])

focs[,5] <- as.character(focs[,5])

#map it with gold standard

#link_frame <- read.table("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/Link_frame/BingRen_Capture_C/0.1_LI_TH")

match_focs <- function(x){
  
  gene_loc <- focs[x,c(1,2,3)]
  
  enh_loc <- focs[x,c(5,6,7)]
  
  ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
              link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))+
  any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
        link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
  
  return(ind)
  
}



match_focs <- Vectorize(match_focs)

foc_ind <- match_focs(1:167988)



sub_info_frame <- subset(info_frame,info_frame$dist<500000&info_frame$act==1)

focs_d <- abs(focs$ep_dist)

my_pair <- list()

#my_pair.2 <- list()


for(i in 1:100){
  #print(i)
  d_s <- (i-1)*5000+1
  
  d_e <- i*5000
  
  count <- sum(focs_d<d_e&focs_d>d_s,na.rm = T)
  
  my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
  
  n <- dim(my_pool)[1]
  
  
  my_pred <- my_pool[order(my_pool$pp,decreasing = F),][1:count,]
  #my_pred.2 <- my_pool[order(my_pool$pp,decreasing = T),][1:count,]
  
  my_pair[[i]] <- my_pred
  #my_pair.2[[i]] <- my_pred.2
  
}

my_pair <- data.table::rbindlist(my_pair)

#my_pair.2 <- data.table::rbindlist(my_pair.2)

my_pair <- as.data.frame(my_pair)

#my_pair.2 <- as.data.frame(my_pair.2)

focs_my_frac <- mean(my_pair[,gd],na.rm = T)

#focs_my_frac.2 <- mean(my_pair.2[,gd],na.rm = T)
}else{
  
  foc_ind <- 0
  
  focs_my_frac <- 1
  
  
}


#compare ernst
print("comparing Ernst")
if(ernst_link!="NA"){
  
  ernst <- read.table(ernst_link)
  
  ernst <- ernst[order(ernst$V4,ernst$V2),]
  
  #merge ernst enhancer
  
  tmp_score<- ernst[1,5]
  
  ernst_merge <- list()
  
  my_merge <- function(ernst){
    
    ernst_merge <- list()
    
    i <- 1
    
    j <- 2
    
    tmp_chr <-  ernst[i,1]
    
    tmp_start <- ernst[i,2]
    
    tmp_stop <- ernst[i,3]
    
    tmp_gene <- ernst[i,4]
    
    tmp_score <- ernst[i,5]
    
    tmp_dist <- ernst[i,6]
    
    nxt_chr <- ernst[j,1]
    
    nxt_start <- ernst[j,2]
    
    nxt_stop <- ernst[j,3]
    
    nxt_gene <- ernst[j,4]
    
    nxt_score <- ernst[j,5]
    
    nxt_dist <- ernst[j,6]
    
    while (!(i > dim(ernst)[1]) & !(j > dim(ernst)[1])) {
      
      
      #if(i%%100==0){print(i)}
      
      while((!tmp_stop<nxt_start)&tmp_chr==nxt_chr&tmp_gene == nxt_gene){
        
        tmp_stop <- nxt_stop
        
        tmp_score <- max(tmp_score,nxt_score)
        
        tmp_dist <- floor(mean(abs(tmp_dist),abs(nxt_dist)))
        
        j <- j+1
        
        nxt_chr <- ernst[j,1]
        
        nxt_start <- ernst[j,2]
        
        nxt_stop <- ernst[j,3]
        
        nxt_gene <- ernst[j,4]
        
        nxt_score <- ernst[j,5]
        
        nxt_dist <- ernst[j,6]
        if(!j<dim(ernst)[1]){break}
        
      }
      
      tmp_frame <- data.frame(tmp_chr,tmp_start,tmp_stop,tmp_gene,tmp_score,tmp_dist)
      
      ernst_merge[[length(ernst_merge)+1]] <- tmp_frame
      
      i <- j
      
      tmp_chr <-  ernst[i,1]
      
      tmp_start <- ernst[i,2]
      
      tmp_stop <- ernst[i,3]
      
      tmp_gene <- ernst[i,4]
      
      tmp_score <- ernst[i,5]
      
      tmp_dist <- ernst[i,6]
      
      j <- j+1
      
      if(!j<dim(ernst)[1]){break}
      nxt_chr <- ernst[j,1]
      
      nxt_start <- ernst[j,2]
      
      nxt_stop <- ernst[j,3]
      
      nxt_gene <- ernst[j,4]
      
      nxt_score <- ernst[j,5]
      
      nxt_dist <- ernst[j,6]
      
    }
    
    ernst_merge <- data.table::rbindlist(ernst_merge)
    
    ernst_merge <- as.data.frame(ernst_merge)
    
    return(ernst_merge)
    
  }
  
  
  ernst_merge <- my_merge(ernst)
  
  #resize ernst enhancer and promoter size
  
  #e_center <- floor((ernst_merge[,2] + ernst_merge[,3])/2)
  
  for(i in 1:dim(ernst_merge)[1]){
    
    
    ernst_merge[i,2] <- ernst_merge[i,2]+50
    ernst_merge[i,3] <- ernst_merge[i,3]-50
    
    
  }
  
  
  #ernst_merge <- my_merge(ernst_merge)
  
  
  # for(i in 1:dim(ernst_merge)[1]){
  #   
  #   
  #   ernst_merge[i,2] <- ernst_merge[i,2]-100
  #   ernst_merge[i,3] <- ernst_merge[i,3]+100
  #   
  #   
  # }
  
  
  
  l <- split(ernst_merge,ernst_merge$tmp_gene)
  
  l <- do.call(rbind,lapply(l, function(x){x[1,]}))
  
  #ernst_merge <- l
  
  ernst_merge <- na.omit(ernst_merge)
  
  gene_annotation <- read.table("../gene_annotation/gene_annotation_V19.txt")
  
  #promoter_frame[,1] <- paste0("chr",promoter_frame[,1])
  
  ernst_merge[,1] <- as.character(ernst_merge[,1])
  
  ernst_merge[,4] <- as.character(ernst_merge[,4])
  
  #gene_annotation[which(gene_annotation[,2]=="23"),2] <- "X"
  
  #tss_loc <- list()
  
  gene_annotation[,7] <- sapply(strsplit(as.character(gene_annotation[,7]),"[.]"),function(x){x[[1]][1]})
  
  
  #for(x in 1:dim(ernst_merge)[1]){
    
    #gene_loc <- gene_annotation[which(gene_annotation[,7]==ernst_merge[x,4]),]
    
    #if(gene_loc[4]==1){
      
      #gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,2]-1000,stop=gene_loc[,2]+1000)
      
    #}else{
      
      #gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
      
    #}
    
    #tss_loc[[i]] <- gene_loc
    
  #}
  
  ernst_merge[,1] <- as.character(ernst_merge[,1])
  
  ernst_merge[,4] <- as.character(ernst_merge[,4])
  
  match_ernst <- function(x){
    
    gene_loc <- gene_annotation[which(gene_annotation[,7]==ernst_merge[x,4]),]
    
    if(dim(gene_loc)[1]==0){return(0)}
    
    if(gene_loc[4]==1){
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,2]-1000,stop=gene_loc[,2]+1000)
      
    }else{
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
      
    }
    
    enh_loc <- ernst_merge[x,c(1,2,3)]
    
    ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
                link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))+p_w*any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
          link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
    
    return(ind)
    
  }
  
  match_ernst <- Vectorize(match_ernst)
  
  ernst_ind <- match_ernst(1:dim(ernst_merge)[1])
  tss <- list()
  
  
  for(i in 1:dim(ernst_merge)[1]){
    
    gene_loc <- gene_annotation[which(gene_annotation[,7]==ernst_merge[i,4]),][1,]
    
    if(any(is.na(gene_loc))){tss[[i]] <- data.frame(chr="NA",start="NA",stop="NA");next}
    
    if(gene_loc[4]==1){
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,2]-1000,stop=gene_loc[,2]+1000)
      
    }else{
      
      gene_loc <- data.frame(chr=as.character(gene_loc[,1]),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
      
    }
    
    tss[[i]] <- gene_loc
    
    
  }
  
  
  tss <- do.call(rbind,tss)
  
  ernst_merge <- cbind(ernst_merge,tss)
  
  
  #sample our pair based on ernst distance
  my_pair <- list()
  
  #my_pair.2 <- list()
  
  ernst_d <- abs(as.numeric(ernst_merge[,8])-as.numeric(ernst_merge[,2]))
  
  sub_info_frame <- subset(info_frame,info_frame$dist<250000&info_frame$act==1)
  
  for(i in 1:100){
    #print(i)
    d_s <- (i-1)*2500+1
    
    d_e <- i*2500
    
    count <- sum(ernst_d<d_e&ernst_d>d_s,na.rm=T)
    
    my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
    my_pred <- my_pool[order(my_pool$corrected_p,decreasing = F),][1:count,]
    #my_pred.2 <- my_pool[order(my_pool$pp,decreasing = T),][1:count,]
    
    my_pair[[i]] <- my_pred
    #my_pair.2[[i]] <- my_pred.2
    
  }
  
  my_pair <- data.table::rbindlist(my_pair)
  
  #my_pair.2 <- data.table::rbindlist(my_pair.2)
  
  my_pair <- as.data.frame(my_pair)
  
  #my_pair.2 <- as.data.frame(my_pair.2)
  
  ernst_my_frac <- mean(my_pair[,gd],na.rm = T)
  
  #ernst_my_frac.2 <- mean(my_pair.2[,gd],na.rm = T)
  
}else{
  
  ernst_ind <- 0
  
  ernst_my_frac <- 1
  
  #ernst_my_frac.2 <- 1
  
}


#target finder
print("comparing TargetFinder")

if(targetfinder_link!="NA"){
  
  targetfinder <- read.table(targetfinder_link,sep=",",header=T)
  
  
  pattern_id <- which(strsplit(as.character(targetfinder[1,1]),"")[[1]]=="|")

  pattern <- substr(targetfinder[1,1],1,pattern_id-1)

  targetfinder_enh <- do.call(rbind,strsplit(gsub(paste0(pattern,"\\|"),"",targetfinder$enhancer_name),"'|'|[:]|-"))
  
  targetfinder_gene <- do.call(rbind,strsplit(gsub(paste0(pattern,"\\|"),"",targetfinder$promoter_name),"'|'|[:]|-"))
  
  targetfinder <- data.frame(targetfinder_enh,targetfinder_gene)
  
  targetfinder[,c(2,3,5,6)] <- t(apply(targetfinder[,c(2,3,5,6)],1,as.numeric))
  
  #resize targetfinder's enhancer to 1000bps
  
  e_center <- floor((targetfinder[,2]+targetfinder[,3])/2)
  
  p_center <- floor((targetfinder[,5]+targetfinder[,6])/2)
  
  targetfinder[,2] <- targetfinder[,2]-400
  
  targetfinder[,3] <- targetfinder[,3]+400
  
  targetfinder[,5] <- p_center-1000
  
  targetfinder[,6] <- p_center+1000
  
  
  match_targetfinder <- function(x){
    
    #gene_loc <- gene_annotation[which(gene_annotation[,1]==targetfinder[x,4]),]
    
    
    
    # if(gene_loc[5]==1){
    #   
    #   gene_loc <- data.frame(chr=as.character(paste0("chr",gene_loc[,2])),start=gene_loc[,3]-1000,stop=gene_loc[,3]+1000)
    #   
    # }else{
    #   
    #   gene_loc <- data.frame(chr=as.character(paste0("chr",gene_loc[,2])),start=gene_loc[,4]-1000,stop=gene_loc[,4]+1000)
    #   
    # }
    
    gene_loc <- targetfinder[x,c(4,5,6)]
    
    enh_loc <- targetfinder[x,c(1,2,3)]
    
    
    ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
                link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))
    any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
          link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
    
    return(ind)
    
  }
  
  match_targetfinder <- Vectorize(match_targetfinder)
  
  targetfinder_ind <- match_targetfinder(1:dim(targetfinder)[1])
  
  
  targetfinder_d <- abs(targetfinder[,5]-targetfinder[,2])
  
  sub_info_frame <- subset(info_frame,info_frame$act==1)
  
  my_pair <- list()
  
  #my_pair.2 <- list()
  
  for(i in 1:100){
    #print(i)
    d_s <- (i-1)*2e4
    
    d_e <- i*2e4
    
    count <- sum(targetfinder_d<d_e&targetfinder_d>d_s)
    
    my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
    my_pred <- my_pool[order(my_pool$corrected_p,decreasing = F),][1:count,]
    #my_pred.2 <- my_pool[order(my_pool$pp,decreasing = T),][1:count,]
    
    my_pair[[i]] <- my_pred
    #my_pair.2[[i]] <- my_pred.2
    
  }
  
  my_pair <- data.table::rbindlist(my_pair)
  
  #my_pair.2 <- data.table::rbindlist(my_pair.2)
  
  my_pair <- as.data.frame(my_pair)
  
  #my_pair.2 <- as.data.frame(my_pair.2)
  
  targetfinder_my_frac <- mean(my_pair[,gd],na.rm = T)
  
  #targetfinder_my_frac.2 <- mean(my_pair.2[,gd],na.rm = T)
  
}else{
  
  targetfinder_ind <- 0
  
  targetfinder_my_frac <- 1
  
  
}

#compare JEME
print("comparing JEME")
if(jeme_link!="NA"){
  
  jeme <- read.csv(jeme_link,header = F)
  
  
  jeme_enh <- as.data.frame(do.call(rbind,strsplit(as.character(jeme[,1]),"[:]|-")))
  
  jeme_gene <- as.data.frame(do.call(rbind,strsplit(as.character(jeme[,2]),"[$]")))
  
  
  jeme <- as.data.frame(cbind(jeme_enh,jeme_gene))
  
  jeme[,c(2,3,7)] <- apply(jeme[,c(2,3,7)],2,as.numeric)
  
  jeme[,c(1,4,5,6)] <- apply(jeme[,c(1,4,5,6)],2,as.character)
  
  
  match_jeme <- function(x){
    
    gene_loc <- jeme[x,c(6,7)]
    
    gene_loc[,2] <- gene_loc[,2]-1000
    
    gene_loc[,3] <- gene_loc[,2]+2000
    
    enh_loc <- jeme[x,c(1,2,3)]
    
    ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
                link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))+
    any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
          link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
    
    return(ind)
    
  }
  
  
  match_jeme <- Vectorize(match_jeme)
  
  jeme_ind <- match_jeme(1:dim(jeme)[1])
  
  jeme_d <- abs(jeme[,7]-jeme[,2])
  
  sub_info_frame <- subset(info_frame,info_frame$act==1)
  
  my_pair <- list()
  
  my_pair.2 <- list()
  
  for(i in 1:100){
    #print(i)
    d_s <- (i-1)*1e4
    
    d_e <- i*1e4
    
    count <- sum(jeme_d<d_e&jeme_d>d_s)
    
    my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
    my_pred <- my_pool[order(my_pool$adj_p,decreasing = F),][1:count,]
    my_pred.2 <- my_pool[order(my_pool$pp,decreasing = T),][1:count,]
    
    my_pair[[i]] <- my_pred
    my_pair.2[[i]] <- my_pred.2
    
  }
  
  my_pair <- data.table::rbindlist(my_pair)
  
  my_pair.2 <- data.table::rbindlist(my_pair.2)
  
  my_pair <- as.data.frame(my_pair)
  
  my_pair.2 <- as.data.frame(my_pair.2)
  
  jeme_my_frac.1 <- mean(my_pair[,gd],na.rm = T)
  
  jeme_my_frac.2 <- mean(my_pair.2[,gd],na.rm = T)
  
}else{
  
  jeme_ind <- 0
  
  jeme_my_frac <- 1
  
}

#compare IM-PET

print("comparing IM_PET")

if(im_pet_link!="NA"){
  
  im_pet <- read.csv(im_pet_link)
  
  im_pet <- im_pet[order(im_pet$Description,im_pet$'IM.PET'),]
  
  im_pet_ec <- (im_pet[,3]+im_pet[,2])/2
  
  im_pet[,2] <- im_pet_ec-1000
  
  im_pet[,3] <- im_pet_ec+1000
  
  match_im_pet <- function(x){
    
    gene_loc <- im_pet[x,c(4,5)]
    
    gene_loc[,2] <- gene_loc[,2]-1000
    
    gene_loc[,3] <- gene_loc[,2]+2000
    
    enh_loc <- im_pet[x,c(1,2,3)]
    
    ind<- any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,4]<as.numeric(enh_loc[,3])&link_frame[,5]>as.numeric(enh_loc[,2])&
                link_frame[,1]==as.character(gene_loc[,1])&link_frame[,2]<as.numeric(gene_loc[,3])&link_frame[,3]>as.numeric(gene_loc[,2]))
    #any(link_frame[,1]==as.character(enh_loc[,1])&link_frame[,2]<as.numeric(enh_loc[,3])&link_frame[,3]>as.numeric(enh_loc[,2])&
          #link_frame[,1]==as.character(gene_loc[,1])&link_frame[,4]<as.numeric(gene_loc[,3])&link_frame[,5]>as.numeric(gene_loc[,2]))
    
    return(ind)
    
  }
  
  
  match_im_pet <- Vectorize(match_im_pet)
  

  id <- which(!duplicated(im_pet[,c(1,2,3,8)]))
  
  im_pet <- im_pet[id,]
  
  im_pet_ind <- match_im_pet(1:dim(im_pet)[1])
  
  im_pet$ind <- im_pet_ind
  
  im_pet_d <- im_pet$DIS
  
  sub_info_frame <- subset(info_frame,info_frame$act==1)
  
  my_pair <- list()
  
  my_pair.2 <- list()
  
  for(i in 1:100){
    #print(i)
    d_s <- (i-1)*2e4
    
    d_e <- i*2e4
    
    count <- sum(im_pet_d<d_e&im_pet_d>d_s)
    
    my_pool <- subset(sub_info_frame,sub_info_frame$dist<d_e&sub_info_frame$dist>d_s)
    
    my_pred <- my_pool[order(my_pool$adj_p,decreasing = F),][1:count,]
    my_pred.2 <- my_pool[order(my_pool$pp,decreasing = T),][1:count,]
    
    my_pair[[i]] <- my_pred
    my_pair.2[[i]] <- my_pred.2
    
  }
  
  my_pair <- data.table::rbindlist(my_pair)
  
  my_pair.2 <- data.table::rbindlist(my_pair.2)
  
  my_pair <- as.data.frame(my_pair)
  
  my_pair.2 <- as.data.frame(my_pair.2)
  
  im_pet_my_frac.1 <- mean(my_pair[,gd],na.rm = T)
  
  im_pet_my_frac.2 <- mean(my_pair.2[,gd],na.rm = T)
  
}else{
  
  im_pet_ind <- 0
  
  im_pet_my_frac <- 1
  
  
  
}




#plot


# #CD4+ Capture C
# mine=1
# 
# ripple_fc <- 0.01572539/0.01723428
# 
# foc_fc <- 0.004024097/0.02727862
# 
# ernst_fc <- 0.06135756/0.0402928
# 
# targetfinder_fc <- 0.03263151/0.03506634
# 
# jeme_fc <- 0.03463491/0.0562336
# 
# im_pet_fc <- 0.02850512/0.03120107
# 
# 
# #CD4+ ChIA-PET
# 
# mine=1
# 
# ripple_fc <- 0.0003671297/0.001029937
# 
# foc_fc <- 5.952806e-05/0.001295353
# 
# ernst_fc <- 0.06135756/0.0402928
# 
# targetfinder_fc <- 0.03263151/0.03506634
# 
# jeme_fc <- 0.03463491/0.0562336
# 
# im_pet_fc <- 0.02850512/0.03120107
# 
# 
# 
# 
# col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")
# 
# fc <- c(mine,ripple_fc,foc_fc,ernst_fc,targetfinder_fc,jeme_fc,im_pet_fc)
# 
# order <- c(1,1+order(fc[-1],decreasing = T))
# 
# name <- c("Integrative\nmodel","Ripple","FOCS","Ernst\net al","TargetFinder","JEME","IM-PET")
# 

#setwd("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/figure/")

col=c("#DB843D","#4473A8","#70588E","#519BAD","#89A54D","#93A8CF","#AA4744")

name <- c("Integrative\nmodel","Ripple","Ernst\net al","FOCS","TargetFinder","JEME","IM-PET")


#adj_p control distance

pth <- "../figures/systematic_conparison_updated/enrichment/"

dir.create(pth,recursive=T)

setwd(pth)

mine=1

ripple_fc <- mean(ripple_ind>0,na.rm=T)/ripple_my_frac

foc_fc <- mean(foc_ind>0,na.rm=T)/focs_my_frac

ernst_fc <- mean(ernst_ind>0,na.rm=T)/ernst_my_frac

targetfinder_fc <- mean(targetfinder_ind>0,na.rm=T)/targetfinder_my_frac

jeme_fc <- mean(jeme_ind>0,na.rm=T)/jeme_my_frac

im_pet_fc <- mean(im_pet_ind>0,na.rm=T)/im_pet_my_frac

#fc <- c(mine,ripple_fc,ernst_fc,foc_fc,targetfinder_fc,jeme_fc,im_pet_fc)

fc <- c(mine,ripple_fc,ernst_fc,foc_fc)
name <- c("Integrative\nmodel","Ripple","Ernst\net al","FOCS")


order <- c(1,1+order(fc[-1],decreasing = T))

pdf(paste0("dist_ctrl_",figure_name,".pdf"))

barplot(fc,col=col,names=name,border=NA,ylab=paste("relative",figure_name,"enrichemnt"),cex.names=0.75)

legend("topright",legend=gsub("\n"," " ,name),fill=col,border = NA,box.lty=0,cex=1.2)

dev.off()


#pp control distance

#setwd("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/figure/systematic_comparison/pp_distance_control/")

#mine=1

#ripple_fc <- mean(ripple_ind>0)/ripple_my_frac.2

#foc_fc <- mean(foc_ind>0)/focs_my_frac.2

#ernst_fc <- mean(ernst_ind>0)/ernst_my_frac.2

#targetfinder_fc <- mean(targetfinder_ind>0)/targetfinder_my_frac.2

#jeme_fc <- mean(jeme_ind)>0/jeme_my_frac.2

#im_pet_fc <- mean(im_pet_ind>0)/im_pet_my_frac.2

#fc <- c(mine,ripple_fc,foc_fc,ernst_fc,targetfinder_fc,jeme_fc,im_pet_fc)

#order <- c(1,1+order(fc[-1],decreasing = T))

#pdf(paste0("post_p_control_distance_comparison_",figure_name,".pdf"))

#barplot(fc[order],col=col[order],names=name[order],border=NA,ylab=paste("relative",figure_name,"enrichemnt"),cex.names=0.75)

#legend("topright",legend=gsub("\n"," " ,name),fill=col,border = NA,box.lty=0,cex=1.2)

#dev.off()



##########################
#all_my_frac.1 <- mean(subset(as.data.frame(info_frame)[,gd],info_frame$pp>0.96))

#all_my_frac.2 <- mean(subset(as.data.frame(info_frame)[,gd],info_frame$adj_p<0.001))


#adj p binarized  control distance

#setwd("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/figure/systematic_comparison/adj_p_binarized/")

#mine=1

#ripple_fc <- mean(ripple_ind>0)/all_my_frac.2

#foc_fc <- mean(foc_ind>0)/all_my_frac.2

#ernst_fc <- mean(ernst_ind>0)/all_my_frac.2

#targetfinder_fc <- mean(targetfinder_ind>0)/all_my_frac.2

#jeme_fc <- mean(jeme_ind>0)/all_my_frac.2

#im_pet_fc <- mean(im_pet_ind>0)/all_my_frac.2

#fc <- c(mine,ripple_fc,foc_fc,ernst_fc,targetfinder_fc,jeme_fc,im_pet_fc)

#order <- c(1,1+order(fc[-1],decreasing = T))


#pdf(paste0("adj_p_binarized_comparison_",figure_name,".pdf"))

#barplot(fc[order],col=col[order],names=name[order],border=NA,ylab=paste("relative",figure_name,"enrichemnt"),cex.names=0.75)

#legend("topright",legend=gsub("\n"," " ,name),fill=col,border = NA,box.lty=0,cex=1.2)

#dev.off()

#pp binarized  control distance

#setwd("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/data/filter_motif/figure/systematic_comparison/pp_binarized/")

#mine=1

#ripple_fc <- mean(ripple_ind>0)/all_my_frac.1

#foc_fc <- mean(foc_ind)>0/all_my_frac.1

#ernst_fc <- mean(ernst_ind>0)/all_my_frac.1

#targetfinder_fc <- mean(targetfinder_ind>0)/all_my_frac.1

#jeme_fc <- mean(jeme_ind>0)/all_my_frac.1

#im_pet_fc <- mean(im_pet_ind>0)/all_my_frac.1

#fc <- c(mine,ripple_fc,foc_fc,ernst_fc,targetfinder_fc,jeme_fc,im_pet_fc)

#order <- c(1,1+order(fc[-1],decreasing = T))

#pdf(paste0("post_p_binarized_comparison_",figure_name,".pdf"))

#barplot(fc[order],col=col[order],names=name[order],border=NA,ylab=paste("relative",figure_name,"enrichemnt"),cex.names=0.75)

#legend("topright",legend=gsub("\n"," " ,name),fill=col,border = NA,box.lty=0,cex=1.1)

#dev.off()
