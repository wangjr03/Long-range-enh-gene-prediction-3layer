#two arguments
#1: enhancer coordinates
#2: chromosome id


args <- commandArgs(TRUE)

################################### load the re-arranged data #############################################

library(data.table)
enh <- fread(args[1])
enh <- as.data.frame(enh)
chrname <- unique(enh$V1)
length <- length(chrname)

filename_mot <- paste0("/mnt/ufs18/rs-027/compbio/wanglab/haowang/Proj6_3_layer_net/hongjie/data/motif_sub/motif_",chrname, ".bed")

i <- as.numeric(args[2]) # i-th chromosome

chrname <- chrname[i]

## load the enhancer and motif in i-th chromosome
enhancer <- subset(enh,enh[,1]==chrname)
motif <- fread(filename_mot[i])

motif <- as.data.frame(motif)

############################# use apply function for every enhancer #####################################

find_motif <- function(one_enh,mot) {
  lower <- as.numeric(one_enh[2])
  upper <- as.numeric(one_enh[3])
  mot_pos <- which( mot$V3>=lower & mot$V4<=upper )
  
  if (length(mot_pos) != 0) {
    mot_inside <- mot[mot_pos,]
    names <- as.matrix(apply(mot_inside[,-5], 1, paste, collapse="_"))
    onename <- apply(names,2,paste, collapse=",")
    return(onename)
  } else {
    return(NA)
  }
  
}

overlap <- apply(enhancer, 1, find_motif, mot=motif)
save_data <- data.frame(gsub(overlap, pattern = " ", replacement = ""))

dir.create("../motif_enh_overlapping")

filename_over <- paste0("../motif_enh_overlapping/overlap_",chrname, ".bed",sep="")
write.table(save_data, file = filename_over, row.names = FALSE, col.names = FALSE, quote = FALSE)

