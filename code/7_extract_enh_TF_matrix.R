#get enh-TF profile, not purified by cell-type specific expression
#1 parameters: chromosome id

args <- commandArgs(TRUE)

i <- as.numeric(args[1]) # i-th chromosome

####################### load the data in i-th chromosome ###############################
library(data.table)
enh <- read.table("../consensus_enhancer_coord.txt")

chrname <- as.character(unique(enh$V1))
length <- length(chrname)

filename_gene <- paste("/mnt/research/compbio/wanglab/hongjie/data/enh_gene_over/overlap_",chrname, ".bed",sep="")
filename_tf <- paste("../motif_enh_overlapping/overlap_",chrname, ".bed",sep="")

## load the enh-gene and enh-TF in i-th chromosome
tf <- read.table(filename_tf[i], sep = "\t",colClass="character")

######################## get enh-TF matrix for each gene ##############################
use_pos <- which(tf != "NA")

tf <- tf[use_pos,]


num_enh <- length(tf)

## define function
get_tf <- function(one_enh) {
  aft_spl <- unlist(strsplit(as.character(one_enh), ",|_"))
  
  name_pos <- seq(1,length(aft_spl),8)
  tf_name <- as.vector(aft_spl[name_pos])
  return(tf_name)
}
get_tf <- Vectorize(get_tf)


## get all TF name list
load("/mnt/research/compbio/wanglab/haowang/Proj6_3_layer_net/preprocess_data/all_tf_name")
tf_name <- unique(all_tf_name)

# #purify un expressed tf
# tf_exp <- read.table("/mnt/research/compbio/wanglab/data/Roadmap/genes/RPKM_all_gene_56epigenomes_simple_select_2.bed")
# 
# ensg_name <- read.table("/mnt/research/compbio/wanglab/data/Roadmap/genes/RPKM_all_gene_name_select_2.bed")
# 
# gene_annotation <- read.table("/mnt/research/compbio/wanglab/data/Roadmap/genes/Gene_annotation.bed")
# 
# tf_ensg_id <- as.character(gene_annotation[match(as.character(tf_name),as.character(gene_annotation[,7])),1])
# 
# tf_loc <- match(tf_ensg_id,as.character(ensg_name[,1]))
# 
# expressed_tf <- subset(tf_name,floor(tf_exp[tf_loc,ct])>0)
# 
# tf_name <- expressed_tf
#

num_tf <- length(tf_name)

## get the enh-TF matrix
enh_tf <- matrix(data = 0, nrow = num_enh, ncol = num_tf)

for (m in 1:num_enh) {
  tf_in_enh <- as.character(unlist(get_tf(tf[m])))
  
  howmany <- length(tf_in_enh)
  
  for (n in 1:howmany) {
    pos <- which(tf_name == tf_in_enh[n])
    enh_tf[m,pos] <- enh_tf[m,pos] + 1
  }
  
}

colnames(enh_tf) <- tf_name


######################## save data for later use ##############################
dir.create("../enh_tf_matrix")
dir.create("../enh_tf_Rdata")

filename_matrix <- paste0("../enh_tf_matrix/enh_tf_",chrname, ".bed",sep="")
write.table(enh_tf, file = filename_matrix[i], row.names = FALSE, col.names = TRUE, quote = FALSE)

filename_2link <- paste0("../enh_tf_Rdata/",chrname, ".Rdata",sep="")

save(
  enh_tf,
  use_pos,
  file = filename_2link[i] 
)





