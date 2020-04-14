args <- commandArgs(T)

ct <- args[1]

gs_name <- args[2]

file.1 <- paste0("../output_promoter/",ct,"/info_frame.Rdata")

load(file.1)

promoter_frame <- info_frame

file.2 <- paste0("../output_exons/",ct,"/info_frame.Rdata")

load(file.2)

exons_frame <- info_frame


promoter_frame <- subset(promoter_frame,promoter_frame$act==1)
exons_frame <- subset(exons_frame,exons_frame$act==1)

library(PRROC)

plot_pr <- function(x,main,col){
  
  plot(x$curve[,c(1,2)],type="l",xlim=c(0,1),ylim=c(0,y_lim),col=col,xlab = "Recall",ylab="Precision")
  
  par(new=T)
  
}



plot_roc <- function(x,main,col){
  
  plot(x$curve[,c(1,2)],type="l",xlim=c(0,1),ylim=c(0,1),col=col,xlab = "FPR",ylab="TPR")
  
  par(new=T)
  
}

gs_promoter <- which(names(promoter_frame)==gs_name)
gs_exons <- which(names(exons_frame)==gs_name)


promoter_roc <- roc.curve(scores.class0 = promoter_frame$corrected_p,weights.class0 = promoter_frame[,gs_promoter],curve=T)
promoter_pr <- pr.curve(scores.class0 = promoter_frame$corrected_p,weights.class0 = promoter_frame[,gs_promoter],curve=T)


exons_roc <- roc.curve(scores.class0 = exons_frame$corrected_p,weights.class0 = exons_frame[,gs_exons],curve=T)
exons_pr <- pr.curve(scores.class0 = exons_frame$corrected_p,weights.class0 = exons_frame[,gs_exons],curve=T)

y_lim <- max(exons_pr$curve[,2],promoter_pr$curve[,2])

dir.create("../RNA_seq_version_comparison/")

setwd("../RNA_seq_version_comparison/")


pdf(paste0(ct,"_",gs_name,".pdf"))

plot_roc(promoter_roc,"","red")
plot_roc(exons_roc,"","orange")

legend("bottomright",legend=c(paste0("Promoter RNA-seq,AUC:",round(promoter_roc$auc,2)),paste0("Exons RNA-seq,AUC:",round(exons_roc$auc,2))),col=c("red","orange"),lty=1)

plot_pr(promoter_pr,"","red")
plot_pr(exons_pr,"","orange")

legend("topright",legend=c(paste0("Promoter RNA-seq,AUPR:",round(promoter_pr$auc.integral,5)),paste0("Exons RNA-seq,AUPR:",round(exons_pr$auc.integral,5))),col=c("red","orange"),lty=1)

dev.off()
