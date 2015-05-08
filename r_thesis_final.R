require(dplyr)
require(tidyr)
require(ggplot2)
require(gridExtra)
require(cowplot)

## Takes a vector of RPKM log2 values (no NAs in nonlogarithmic data)
countNAs <- function(x){
  return (length(x) - length(na.exclude(x)))
}

isComRPKM <- function(chi){
  return (!(0.0 %in% chi))
}

is_rpkm_complete <- function(x){
  return (!(0.0 %in% x))
}
## split (HE/LE KDE)
is_complete <- function(x){
  return (countNAs(x) == 0)
}

is_high <- function(x, lim=5.0){
  return (!(mean(x, na.rm=T) <= lim))
}

remove_non_tumor <- function(x){
  if (length(x %>% select(contains("HYB"))) == 1){
    colnames(x)[1] <- "gene"
  }
  gene_names <- x %>% select(gene)
  vals <- x %>% select(contains(".01"))
  ret <- cbind(gene_names, vals)
  return (as_data_frame(ret))
}

## Takes a transposed, log2 RPKM file
make_numNAs_on_means_plots <- function(chi){
chi_no_gene <- chi %>% select(-contains("gene"))

x <- apply(chi_no_gene, 2, mean, na.rm=T)
n <- colnames(chi_no_gene)
y <- apply(chi_no_gene, 2, countNAs)

mac <- as_data_frame(data.frame(nam = n, means = x, numNAs = y))
mac %>% arrange(desc(means)) -> mac
return (qplot(mac$means, mac$numNAs, xlim = c(-5, 20)), xlab = "Mean of miRNA Expression")
#print(qplot(mac$means, mac$numNAs, xlim = c(-5, 20), ylim=c(0,5)) "Number of NA Values")
}

make_tri_plot <- function(chi, title, s_ylim=c(0,800)){
  x <- chi %>% select(-gene)
  x_full <- unlist(x)
  p1 <- qplot(x_full, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab="All Genes (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  x_missing <- unlist(x[!complete.cases(x),])
  p2 <- qplot(x_missing, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab="Genes Missing in >= 1 Sample (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  x_no_missing <- unlist(x[complete.cases(x),])
  p3 <- qplot(x_no_missing, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab = "Genes Present in All Samples (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  plts <- list(p1, p2, p3)
  return (plot_grid(plts[[1]], plts[[2]], plts[[3]], rows=3, labels=c("A", "B", "C")))
}

make_L_plots_dplyr <- function(chi, title){
  par(ask=T)
  chi <- chi %>% gather(gene, "value")
  colnames(chi)[2] <- "sample"
  chi <- chi %>% group_by(gene) %>% mutate(numNAs=countNAs(value))
  chi <- chi %>% group_by(gene) %>% mutate(med=median(value, na.rm=T))
  chi <- chi %>% group_by(gene) %>% mutate(mad=mad(value, na.rm=T))
  print(ggplot(chi, aes(x=med, y=numNAs)) +
          geom_point(aes(color=mad, size=mad), alpha=.9) +
          theme(axis.title.y = element_text(vjust=-0.001)) + xlim(c(-5, 20)) +
          xlab("Number of Samples Missing miRNA") +
          ylab("Median Log2 RPKM of miRNA") + ggtitle(title))
  print(ggplot(chi, aes(x=med, y=numNAs)) +
                geom_point(aes(color=mad, size=mad), alpha=.9) +
                theme(axis.title.y = element_text(vjust=-0.001)) + xlim(c(-5, 20)) + ylim(c(0, 10)) +
                xlab("Median Log2 RPKM of miRNA") + ylab("Number of Samples Missing miRNA") + ggtitle(title)
  )
}


make_qqnorm <- function(x,y, title = "QQ Plot"){
  qqplot(x,y, main = title)
}

fit_distribs_helper <- function(d){
  ## A fitting function borrowed from World of Piggy
  tmp <- as.numeric(d)
  tmp <- tmp[!is.na(tmp)]
  ret <- c()
  for (i in c("normal", "gamma")){
    fit <- fitdistr(tmp, i)
    est_shape <- fit$estimate[[1]]
    est_rate <- fit$estimate[[2]]
    ks = ks.test(tmp, "pgamma", shape=est_shape, rate=est_rate)
    
  }
}

make_cor_plots_median_log2 <- function(chi){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=mad(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing_Samples=is_complete(value))
  return (ggplot(tmp) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing_Samples)) +
    labs(x="Median(Log2 RPKM)", y="MAD(Log2 RPKM)") + xlim(c(0,20)))
}

make_cor_plots_mean_log2 <- function(chi, title){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=mean(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=sd(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing_Samples=is_complete(value))
  ggplot(tmp) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing_Samples)) +
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Mean(RPKM)", y="SD(RPKM)") + xlim(c(0,20)) + ggtitle(title)#+ geom_abline(aes(color="red"))
}

make_cor_plots_untrans_median <- function(chi){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(devi=mad(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(isCom=isComRPKM(value))
  ggplot(tmp) + geom_point(aes(x=avg, y=devi, size=devi, color=isCom)) +
    labs(x="Median(RPKM)", y="MAD(RPKM)") + geom_abline()
  #     scale_size_continuous(name="MAD(RPKM)", breaks=c(25000, 50000, 75000 ), labels=c("25,000", "50,000", "75,000"))
}

make_cor_plots_untrans_mean <- function(chi, title){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=mean(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=sd(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing=isComRPKM(value))
  ggplot(tmp) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing)) +
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Mean(RPKM)", y="SD(RPKM)") + geom_abline() + ggtitle(title)
}

make_cor_plots_untrans_median_zoom <- function(chi){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=mad(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing=isComRPKM(value))
  return (ggplot(tmp %>% group_by(gene)) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing)) +
    labs(x="Median(RPKM)", y="MAD(RPKM)") + geom_abline() + xlim(0, 100) + ylim(c(0, 100)))
}


ind_kde_plot <- function(x, k=.10){
  y <- x %>% sample_frac(k) %>% gather(gene)
  ret <- ggplot(y, aes(x = value, fill = gene)) + geom_density(alpha = .2) + labs(x = "Log2 RPKM") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0)
  return (ret)
}


make_KDE_plot <- function(chi, title, k=100, individual=F){
#   samples <<- chi %>% sample_frac(k, replace = FALSE)
    #zed <<- samples %>% select(-contains("gene"))
    chi <- as.data.frame(chi, row.names = chi$gene)
  if (individual == T){
    y <- data.frame()
    samples <- sample(rownames(chi), k)
    for (i in samples){
      x <- data.frame(vals = t(chi[i,]), row.names=NULL)
      colnames(x)[1] <- "vals"
      x$samp <- i
      y <- rbind(y,x)
    }
    
    return (ggplot(y, aes(vals)) + geom_density(alpha = 0.4)  + ggtitle(title) + labs(x = "Log(2) RPKM") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0))
  }
  else{
    z_tot <<- unlist(zed)
#     print(ggplot(zed, aes(x = zed)) + geom_density(alpha = 0.4) + scale_fill_hue(l=30) + ggtitle(title) + labs(x = "RPKM (Log2") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0))

  }
}

# acc$c_type <- "ACC"
# brca$c_type <- "BRCA"
# combi <- bind_rows(acc, brca)
# ggplot(combi, aes(x=value, fill=c_type)) + geom_density(alpha=.4) + scale_fill_discrete(name="Cancer Type") + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for ACC and BRCA")
# 
# 
# brca_full <- as_data_frame(data.frame(vals=unlist(BRCA.RPKM.log2.munged), sample="BRCA"))
# 
# ggplot(brca, aes(x=value)) + 
#   geom_density(aes(x=value, colour="green"), alpha=0.4) +
#   geom_density(aes(x=value, fill=isCom), alpha=0.2) +
#   xlim(-5,20) + ggtitle("KDE for BRCA") + ylim(0, 1.0)+
# theme(axis.title.y = element_text(vjust=-0.001)) + theme(legend.position="none")
# legend("topright", c("Overall", "miRNAs w/ Missing Samples", "miRNAs In All Samples"), col=c(2:4), pch=19, cex=0.4)
# 
# ggplot(brca, aes(x=value, fill=gene)) + geom_density(alpha=.4)  + xlim(c(-5, 20)) + ylim(c(0,1.0))  + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for BRCA, miRNA") + theme(legend.position="none")



# log_mad <- function(x, na.rm = F){
#   ret <- log(mad(x, na.rm))
#   return(ret)
# }
# x <- apply(ACC.RPKM.log2.munged.transposed[,-1], 2, mean, na.rm=T)
# y <- apply(ACC.RPKM.log2.munged.transposed[,-1], 2, log_mad, na.rm=T)
# qplot(x, y, main = "Log(MAD) on Log(Mean) of Expression")









get_lost_expr <- function(chi){
  chi %>% gather(gene, "value") -> chi
  colnames(chi)[2] <- "sample"
  chi %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T)) -> chi
  chi %>% group_by(gene) %>% mutate(numNAs=countNAs(value)) -> chi
  ret <- chi %>% group_by(gene) %>% filter(avg > 5.0) %>% filter(numNAs > 0)
  return (ret)
}


# for (i in c(1:14)){
#   #   cat(str(i))
#   #   cat("\n")
#   #   colnames(all_log2_frames[i])[1] <- "gene"
# #   message(log_names[i])
#   tmp <- get_lost_expr(all_log2_frames[[i]])
#   ret <- spread(tmp, sample, value)
#   message(log_names[i])
#   message(paste(ret$gene, sep=" ", collapse=" "))
# }


##
## Dendrograms, but with only LE gene




## Raw RPKM
ACC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ACC.RPKM.munged.txt", sep="\t"))
ACC.RPKM.munged <- remove_non_tumor(ACC.RPKM.munged)
acc <- prep_for_dendro(ACC.RPKM.munged, 1)

BLCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BLCA.RPKM.munged.txt", sep="\t"))
BLCA.RPKM.munged <- remove_non_tumor(BLCA.RPKM.munged)
blca <- prep_for_dendro(BLCA.RPKM.munged, 2)

BRCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.munged.txt", sep="\t"))
BRCA.RPKM.munged <- remove_non_tumor(BRCA.RPKM.munged)
brca <- prep_for_dendro(BRCA.RPKM.munged, 3)

CESC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/CESC.RPKM.munged.txt", sep="\t"))
CESC.RPKM.munged <- remove_non_tumor(CESC.RPKM.munged)
cesc <- prep_for_dendro(CESC.RPKM.munged, 4)

COADREAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COADREAD.RPKM.munged.txt", sep="\t"))
COADREAD.RPKM.munged <- remove_non_tumor(COADREAD.RPKM.munged)
coadread <- prep_for_dendro(COADREAD.RPKM.munged, 5)

COAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COAD.RPKM.munged.txt", sep="\t"))
COAD.RPKM.munged <- remove_non_tumor(COAD.RPKM.munged)
coad <- prep_for_dendro(COAD.RPKM.munged, 6)

DLBC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/DLBC.RPKM.munged.txt", sep="\t"))
DLBC.RPKM.munged <- remove_non_tumor(DLBC.RPKM.munged)
dlbc <- prep_for_dendro(DLBC.RPKM.munged, 7)

ESCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ESCA.RPKM.munged.txt", sep="\t"))
ESCA.RPKM.munged <- remove_non_tumor(ESCA.RPKM.munged)
esca <- prep_for_dendro(ESCA.RPKM.munged, 8)

HNSC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/HNSC.RPKM.munged.txt", sep="\t"))
HNSC.RPKM.munged <- remove_non_tumor(HNSC.RPKM.munged)
hnsc <- prep_for_dendro(HNSC.RPKM.munged, 9)

KICH.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KICH.RPKM.munged.txt", sep="\t"))
KICH.RPKM.munged <- remove_non_tumor(KICH.RPKM.munged)
kich <- prep_for_dendro(KICH.RPKM.munged, 10)

KIRC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KIRC.RPKM.munged.txt", sep="\t"))
KIRC.RPKM.munged <- remove_non_tumor(KIRC.RPKM.munged)
kirc <- prep_for_dendro(KIRC.RPKM.munged, 11)

LAML.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/LAML.RPKM.munged.txt", sep="\t"))
laml <- prep_for_dendro(LAML.RPKM.munged, 12)

PAAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PAAD.RPKM.munged.txt", sep="\t"))
PAAD.RPKM.munged <- remove_non_tumor(PAAD.RPKM.munged)
paad <- prep_for_dendro(PAAD.RPKM.munged, 13)

PRAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PRAD.RPKM.munged.txt", sep="\t"))
PRAD.RPKM.munged <- remove_non_tumor(PRAD.RPKM.munged)
prad <- prep_for_dendro(PRAD.RPKM.munged, 14)




## log2 Expression, for tri plots
ACC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ACC.RPKM.log2.munged.txt", sep="\t"))
ACC.RPKM.log2.munged <- remove_non_tumor(ACC.RPKM.log2.munged)


BLCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BLCA.RPKM.log2.munged.txt", sep="\t"))
BLCA.RPKM.log2.munged <- remove_non_tumor(BLCA.RPKM.log2.munged)

BRCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.log2.munged.txt", sep="\t"))
BRCA.RPKM.log2.munged <- remove_non_tumor(BRCA.RPKM.log2.munged)

CESC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/CESC.RPKM.log2.munged.txt", sep="\t"))
CESC.RPKM.log2.munged <- remove_non_tumor(CESC.RPKM.log2.munged)

COADREAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COADREAD.RPKM.log2.munged.txt", sep="\t"))
COADREAD.RPKM.log2.munged <- remove_non_tumor(COADREAD.RPKM.log2.munged)

COAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COAD.RPKM.log2.munged.txt", sep="\t"))
COAD.RPKM.log2.munged <- remove_non_tumor(COAD.RPKM.log2.munged)

DLBC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/DLBC.RPKM.log2.munged.txt", sep="\t"))
DLBC.RPKM.log2.munged <- remove_non_tumor(DLBC.RPKM.log2.munged)

ESCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ESCA.RPKM.log2.munged.txt", sep="\t"))
ESCA.RPKM.log2.munged <- remove_non_tumor(ESCA.RPKM.log2.munged)

HNSC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/HNSC.RPKM.log2.munged.txt", sep="\t"))
HNSC.RPKM.log2.munged <- remove_non_tumor(HNSC.RPKM.log2.munged)

KICH.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KICH.RPKM.log2.munged.txt", sep="\t"))
KICH.RPKM.log2.munged <- remove_non_tumor(KICH.RPKM.log2.munged)

KIRC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KIRC.RPKM.log2.munged.txt", sep="\t")) 
KIRC.RPKM.log2.munged <- remove_non_tumor(KIRC.RPKM.log2.munged)

LAML.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/LAML.RPKM.log2.munged.txt", sep="\t"))

PAAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PAAD.RPKM.log2.munged.txt", sep="\t"))
PAAD.RPKM.log2.munged <- remove_non_tumor(PAAD.RPKM.log2.munged)

PRAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PRAD.RPKM.log2.munged.txt", sep="\t"))
PRAD.RPKM.log2.munged <- remove_non_tumor(PRAD.RPKM.log2.munged)

dframes_log <- list(
  ACC.RPKM.log2.munged, BLCA.RPKM.log2.munged, BRCA.RPKM.log2.munged,
  COAD.RPKM.log2.munged, COADREAD.RPKM.log2.munged, CESC.RPKM.log2.munged, DLBC.RPKM.log2.munged,
  ESCA.RPKM.log2.munged, HNSC.RPKM.log2.munged, KICH.RPKM.log2.munged, KIRC.RPKM.log2.munged, LAML.RPKM.log2.munged,
  PAAD.RPKM.log2.munged, PRAD.RPKM.log2.munged)

dframes_rpkm <- list(
  ACC.RPKM.munged, BLCA.RPKM.munged, BRCA.RPKM.munged,
  COAD.RPKM.munged, COADREAD.RPKM.munged, CESC.RPKM.munged, DLBC.RPKM.munged,
  ESCA.RPKM.munged, HNSC.RPKM.munged, KICH.RPKM.munged, KIRC.RPKM.munged, LAML.RPKM.munged,
  PAAD.RPKM.munged, PRAD.RPKM.munged)


log_names <- list("ACC", "BLCA", "BRCA",
                  "COAD", "COADREAD", "CESC", "DLBC",
                  "ESCA", "HNSC", "KICH", "KIRC", "LAML",
                  "PAAD", "PRAD")

for (i in c(1:14)){
  message(log_names[[i]])
  dframes_log[[i]] <- dframes_log[[i]] %>% mutate(ctype=log_names[[i]])
}

for (i in c(1:14)){
  message(log_names[[i]])
  dframes_rpkm[[i]] <- dframes_rpkm[[i]] %>% mutate(ctype=log_names[[i]])
}
## Figure 1, BRCA Tri Plot
brca <- dframes_log[[3]]
x <- make_tri_plot(brca %>% select(-contains("ctype")), "BRCA", c(0, 800))
x

## Figure 2, PanCan Tri Plot
tmp_fin <- do.call(bind_rows, dframes_log)
tmp_fin <- tmp_fin %>% gather(sample, value, -ctype, -gene)
#TODO
a <- ggplot(tmp_fin) + geom_histogram(aes(x=value), binwidth=.3) + xlim(c(-5, 20)) + xlab("Log2 RPKM")
a + facet_wrap(~ ctype, ncol = 3, scales = "free_y")

## Figure 3, Overall mRNA KDE
BRCA.mRNA <- as_data_frame(read.csv("/home/eric/sync/data/brca/BRCA.mRNAseq_RPKM_log2.txt", sep="\t"))
colnames(BRCA.mRNA)[1] <- "gene"
BRCA.mRNA <- remove_non_tumor(BRCA.mRNA)
BRCA.mRNA <- BRCA.mRNA %>% gather(gene)
colnames(BRCA.mRNA)[2] <- "sample"
BRCA.mRNA = BRCA.mRNA %>% group_by(gene) %>% mutate(isComplete=is_complete(value))
a <- qplot(value, ..count.., data=BRCA.mRNA %>% sample_frac(.50), geom="density",adjust=1, xlim=c(-5, 20), xlab="log2 RPKM")
b <- qplot(value, ..count.., data=BRCA.mRNA %>% sample_frac(.50), geom="density", color=isComplete, adjust=1, xlim=c(-5, 20), xlab="log2 RPKM")
plot_grid(a, b, rows=2, labels=c("A", "B"))

## Figure 4, Hebenstreit Figure 1

## Figure 5, Individual KDE for mRNA
BRCA.mRNA <- as_data_frame(read.csv("/home/eric/sync/data/brca/BRCA.mRNAseq_RPKM_log2.txt", sep="\t"))
colnames(BRCA.mRNA)[1] <- "gene"
BRCA.mRNA <- remove_non_tumor(BRCA.mRNA)
a <- ind_kde_plot(BRCA.mRNA, .05)
a



## Figure 6, Overall microRNA KDE,
## Missing/Nonmissing KDE, Missing nonmissing/missing + LE/HE KDE,
## Composite KDE for microRNA
## All for BRCA
## KDE of Total Dists
brcaMI <- gather(BRCA.RPKM.log2.munged, gene, "value")
colnames(brcaMI)[2] <- "sample"
brcaMI <- brcaMI %>% group_by(gene) %>% mutate(isCom=is_complete(value))
brcaMI <- brcaMI %>% group_by(gene) %>% mutate(isHigh=is_high(value, 4.0))
#a <- ggplot(brca, aes(x=value, fill=isHigh)) + geom_density(alpha=.4)  + xlim(c(-5, 20)) + ylim(c(0,.5)) + scale_fill_discrete(name="AVG Log2(RPKM)", breaks=c("TRUE", "FALSE"), labels=c(">5 RPKM", "<=5 RPKM")) + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for BRCA")
a <- qplot(value, ..count.., data=brcaMI %>% sample_frac(.50), geom="density",adjust=1, xlim=c(-5, 20))
b <- ggplot(brcaMI, aes(x=value, fill=isCom)) + geom_density(alpha=.4)  +
  xlim(c(-5, 20)) + ylim(c(0,.5)) +
  scale_fill_discrete(name="Missing Values", breaks=c("TRUE", "FALSE"), labels=c("None", "At least one sample missing")) #+
  #theme(axis.title.y = element_text(vjust=-0.001))
c <- ggplot(brcaMI) +
  geom_density(alpha=.4, aes(x=value, fill=isHigh)) +
  geom_density(alpha=.4, aes(x=value, color=isCom)) +
  xlim(c(-5, 20)) + ylim(c(0,.5)) # +
  #theme(axis.title.y = element_text(vjust=-0.001))
d <- ind_kde_plot(BRCA.RPKM.log2.munged, 1.0)
plot_grid(a, b, c, d, rows=2, labels=c("A", "B", "C", "D"))


## Figure 7, MAD on Median RPKM, Zoomed MAD on Median RPKM, Log2 MAD on Median
##
## A plot of number of missing values on mean RPKM
## shows that above 4-6RPKM, there are no more missing values.
##
a <- make_cor_plots_untrans_median(BRCA.RPKM.munged)
b <- make_cor_plots_untrans_median_zoom(BRCA.RPKM.munged)
c <- make_cor_plots_median_log2(BRCA.RPKM.log2.munged)
plot_grid(a, b, c, rows=3, labels=c("A", "B", "C"))


## Figure 8, BRCA missing on med
BRCA.RPKM.log2.munged.transposed <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.log2.munged.txt.transposed", sep="\t"))
BRCA.RPKM.log2.munged.transposed_no_gene <- BRCA.RPKM.log2.munged.transposed %>% select(-contains("gene"))

x <- apply(BRCA.RPKM.log2.munged.transposed_no_gene, 2, median, na.rm=T)
n <- colnames(BRCA.RPKM.log2.munged.transposed_no_gene)
y <- apply(BRCA.RPKM.log2.munged.transposed_no_gene, 2, countNAs)

mac <- as_data_frame(data.frame(nam = n, medians = x, numNAs = y))
mac %>% arrange(desc(medians)) -> mac
a <- qplot(mac$medians, mac$numNAs, xlim = c(-5, 20)) +
  xlab("Median Log2 RPKM") + ylab("Number of Samples Missing microRNA")

b <- qplot(mac$medians, mac$numNAs, xlim = c(-5, 20), ylim=c(0,5)) +
  xlab("Median Log2 RPKM") + ylab("Number of Samples Missing miRNA")  + ggtitle("(BRCA)")

plot_grid(a, b, rows=1, labels=c("A", "B"))

## Figure 9, PanCan missing on med
tmp_fin_mm <- tmp_fin %>% group_by(gene) %>% mutate(deviation=mad(value, na.rm=T))
tmp_fin_mm <- tmp_fin_mm %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
tmp_fin_mm <- tmp_fin_mm %>% group_by(gene) %>% mutate(numNAs=countNAs(value))

ttt <- BRCA.RPKM.log2.munged %>% gather(gene, "value")
colnames(ttt)[2] <- "sample"
ttt <- ttt %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
ttt <- ttt %>% group_by(gene) %>% mutate(numNAs=countNAs(value))
cor.test(ttt$value, ttt$numNAs)

a <- ggplot(tmp_fin_mm) + geom_point(aes(x=numNAs, y=avg, size=deviation, color=deviation))
a + facet_wrap(~ ctype, ncol = 3, scales = "free_y")

## Figure 10, PanCan scaled Missing on Med
facet_wrap()
## Figure 11, Dendrograms
##
## Raw Expression dendrograms
## for both LE and HE expressed genes.
par(mfrow=c(3,1))
prep_for_dendro<- function(chi, type_val, check_complete=F){
  colnames(chi)[1] <- "gene"
  chi <- chi %>% gather(gene, "value") 
  colnames(chi)[2] <- "sample"
  chi <- chi %>% group_by(gene) %>% mutate(complete=is_rpkm_complete(value))
  if (check_complete){
    chi <- chi %>% filter(complete==T)
  }
  chi <- chi %>% select(-complete) %>% spread(gene, value)
  chi <- chi %>% mutate(type=type_val)
  return (chi)
}
prep_for_dendroLE<- function(chi, type_val, le=T){
  colnames(chi)[1] <- "gene"
  chi <- chi %>% gather(gene, "value")
  colnames(chi)[2] <- "sample"
  chi <- chi %>% group_by(gene) %>% mutate(complete=is_rpkm_complete(value))
  if (le == T){
    chi <- chi %>% filter(complete == F)
  }
  else{
    chi <- chi %>% group_by(gene) %>% filter(complete == T)
  }
  chi <- chi %>% select(-complete) %>% spread(gene, value)
  chi <- chi %>% mutate(type=type_val)
  return (chi)
}


acc <- prep_for_dendro(ACC.RPKM.munged, 1)
blca <- prep_for_dendro(BLCA.RPKM.munged, 2)
brca <- prep_for_dendro(BRCA.RPKM.munged, 3)
cesc <- prep_for_dendro(CESC.RPKM.munged, 4)
coadread <- prep_for_dendro(COADREAD.RPKM.munged, 5)
coad <- prep_for_dendro(COAD.RPKM.munged, 6)
dlbc <- prep_for_dendro(DLBC.RPKM.munged, 7)
esca <- prep_for_dendro(ESCA.RPKM.munged, 8)
hnsc <- prep_for_dendro(HNSC.RPKM.munged, 9)
kich <- prep_for_dendro(KICH.RPKM.munged, 10)
kirc <- prep_for_dendro(KIRC.RPKM.munged, 11)
laml <- prep_for_dendro(LAML.RPKM.munged, 12)
paad <- prep_for_dendro(PAAD.RPKM.munged, 13)
prad <- prep_for_dendro(PRAD.RPKM.munged, 14)

require(sparcl)
tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))

clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab = "", sub = "")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=.8)

##
## Dendros for just LE genes
##

acc <- prep_for_dendroLE(ACC.RPKM.munged, 1)
blca <- prep_for_dendroLE(BLCA.RPKM.munged, 2)
brca <- prep_for_dendroLE(BRCA.RPKM.munged, 3)
cesc <- prep_for_dendroLE(CESC.RPKM.munged, 4)
coadread <- prep_for_dendroLE(COADREAD.RPKM.munged, 5)
coad <- prep_for_dendroLE(COAD.RPKM.munged, 6)
dlbc <- prep_for_dendroLE(DLBC.RPKM.munged, 7)
esca <- prep_for_dendroLE(ESCA.RPKM.munged, 8)
hnsc <- prep_for_dendroLE(HNSC.RPKM.munged, 9)
kich <- prep_for_dendroLE(KICH.RPKM.munged, 10)
kirc <- prep_for_dendroLE(KIRC.RPKM.munged, 11)
laml <- prep_for_dendroLE(LAML.RPKM.munged, 12)
paad <- prep_for_dendroLE(PAAD.RPKM.munged, 13)
prad <- prep_for_dendroLE(PRAD.RPKM.munged, 14)


tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))
clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab = "", sub = "")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.8)


##
## Dendros for HE
##
acc <- prep_for_dendroLE(ACC.RPKM.munged, 1, F)
blca <- prep_for_dendroLE(BLCA.RPKM.munged, 2, F)
brca <- prep_for_dendroLE(BRCA.RPKM.munged, 3, F)
cesc <- prep_for_dendroLE(CESC.RPKM.munged, 4, F)
coadread <- prep_for_dendroLE(COADREAD.RPKM.munged, 5, F)
coad <- prep_for_dendroLE(COAD.RPKM.munged, 6, F)
dlbc <- prep_for_dendroLE(DLBC.RPKM.munged, 7, F)
esca <- prep_for_dendroLE(ESCA.RPKM.munged, 8, F)
hnsc <- prep_for_dendroLE(HNSC.RPKM.munged, 9, F)
kich <- prep_for_dendroLE(KICH.RPKM.munged, 10, F)
kirc <- prep_for_dendroLE(KIRC.RPKM.munged, 11, F)
laml <- prep_for_dendroLE(LAML.RPKM.munged, 12, F)
paad <- prep_for_dendroLE(PAAD.RPKM.munged, 13, F)
prad <- prep_for_dendroLE(PRAD.RPKM.munged, 14, F)

tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))
tmp_fin[is.na(tmp_fin)] <- 0.0
clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab = "", sub = "")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.6)


##
## Distribution fitting
##

## Overall RPKM
brca <- BRCA.RPKM.munged %>% gather(gene, value)
colnames(brca)[2] <- "sample"
qqnorm(brca$value)
abline(0, 1)

## HE RPKM
brcaHE <- brca %>% group_by(gene) %>% mutate(isCom=isComRPKM(value))
brcaHE <- brcaHE %>% group_by(gene) %>% filter(isCom == T)
qqnorm(brcaHE$value)
abline(0, 1)

## LE RPKM
brcaLE <- brca %>% group_by(gene) %>% mutate(isCom=isComRPKM(value))
brcaLE <- brcaLE %>% group_by(gene) %>% filter(isCom == F)
qqnorm(brcaLE$value)
abline(0, 1)

## Overall Log2
brca_log <- BRCA.RPKM.log2.munged %>% gather(gene, value)
colnames(brca_log)[2] <- "sample"
qqnorm(brca_log$value)
abline(0,1)

brca_logHE <- brca_log %>% group_by(gene) %>% mutate(isCom=is_complete(value))
brca_logHE <- brca_logHE %>% group_by(gene) %>% filter(isCom == T)
qqnorm(brca_logHE$value)
abline(0, 1)

brca_logLE <- brca_log %>% group_by(gene) %>% mutate(isCom=is_complete(value))
brca_logLE <- brca_logLE %>% group_by(gene) %>% filter(isCom == F)
qqnorm(brca_logLE$value)
abline(0, 1)

## Set up distributions for fitting
require(dplyr)
require(MASS)
brca_log
x_logHE <- as.numeric(brca_logHE$value)
x_logLE <- na.omit(as.numeric(brca_logLE$value))
fitLE <- fitdistr(x_logLE, "normal")
fitHE <- fitdistr(x_logHE, "normal")
fitLE$estimate
fitLE$n
fitHE$estimate
fitHE$n
ks.test(rnorm(fitHE$n, fitHE$estimate[1], fitHE$estimate[2]), x_logHE)

x_rpkmHE <- as.numeric(brcaHE$value)


