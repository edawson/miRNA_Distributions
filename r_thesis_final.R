require(dplyr)
require(tidyr)
require(ggplot2)
require(gridExtra)
require(cowplot)

## Takes a vector of RPKM log2 values (no NAs in nonlogarithmic data)
countNAs <- function(x){
  return (length(x) - length(na.exclude(x)))
}

remove_non_tumor <- function(x){
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
print(qplot(mac$means, mac$numNAs, xlim = c(-5, 20)), xlab = "Mean of miRNA Expression")
print(qplot(mac$means, mac$numNAs, xlim = c(-5, 20), ylim=c(0,5)) "Number of NA Values")
}


make_tri_plot <- function(chi, title, s_ylim=c(0,800)){
  x <- chi %>% select(-gene)
  x_full <- unlist(x)
  p1 <- qplot(x_full, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab="All Genes (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  x_missing <- unlist(x[!complete.cases(x),])
  p2 <- qplot(x_missing, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab="Genes Missing in >= 1 Sample (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  x_no_missing <- unlist(x[complete.cases(x),])
  p3 <- qplot(x_no_missing, geom = "histogram", binwidth=.01, xlim=c(-5,20), ylim=s_ylim, xlab = "Genes Present in All Samples (Log2 RPKM)", ylab="Count") + theme(axis.title.y = element_text(vjust=-0.001))
  
  ret <- grid.arrange(p1, p2, p3, main = title)
#   return (x)
  return (ret)
}

make_L_plots_dplyr <- function(chi, title){
  par(ask=T)
  chi <- chi %>% gather(gene, "value")
  colnames(chi)[2] <- "sample"
  chi <- chi %>% group_by(gene) %>% mutate(numNAs=countNAs(value))
  chi <- chi %>% group_by(gene) %>% mutate(med=median(value, na.rm=T))
  chi <- chi %>% group_by(gene) %>% mutate(mad=mad(value, na.rm=T))
#   print(ggplot(chi, aes(x=numNAs, y=med)) +
#           geom_point(aes(color=mad, size=mad), alpha=.9) +
#           theme(axis.title.y = element_text(vjust=-0.001)) + ylim(c(-5, 20)) +
#           xlab("Number of Samples Missing miRNA") +
#           ylab("Median Log2 RPKM of miRNA")
#   )
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
#   ggplot(x_noNAs, aes(x=x_noNAs$numberNAs, y=x_noNAs$mean)) +
#     geom_point(aes(color=factor(x_noNAs$nonas), size=x_noNAs$mad), alpha=.5) +
#     theme(axis.title.y = element_text(vjust=-0.001)) + ylim(c(-5, 20))+
#     xlab("Number of Samples Missing miRNA") +
#     ylab("Mean Log2 RPKM of miRNA")
}



make_L_plots <- function(x, title){
  par(ask=T)
  x$numberNAs <- apply(x, 1, countNAs)
  x$numberValues <- apply(select(x, contains("TCGA")), 1, length)
  x$mean <- apply(select(x, contains("TCGA")), 1, mean, na.rm = T)
  x$median <- apply(select(x, contains("TCGA")), 1, median, na.rm = T)
  x$mad <- apply(select(x, contains("TCGA")), 1, mad, na.rm = T)
  x_noNAs <<- x %>% filter(mean != "NA", mad != "NA", median != "NA")
  
  x_noNAs <<- x_noNAs %>% filter(numberNAs / numberValues < .99)
  cor.test(x_noNAs$numberNAs, x_noNAs$mean)
  x_noNAs$nonas <<- x_noNAs$numberNAs == 0
  print(ggplot(x_noNAs, aes(x=x_noNAs$numberNAs, y=x_noNAs$median)) +
    geom_point(aes(color=x_noNAs$mad, size=x_noNAs$mad), alpha=.9) +
    theme(axis.title.y = element_text(vjust=-0.001)) + ylim(c(-5, 20)) +
    xlab("Number of Samples Missing miRNA") +
    ylab("Median Log2 RPKM of miRNA") + ggtitle(title) + scale_y_discrete()
  )
  ggplot(x_noNAs, aes(x=x_noNAs$numberNAs, y=x_noNAs$mean)) +
    geom_point(aes(color=factor(x_noNAs$nonas), size=x_noNAs$mad), alpha=.5) +
    theme(axis.title.y = element_text(vjust=-0.001)) + ylim(c(-5, 20))+
    xlab("Number of Samples Missing miRNA") +
    ylab("Median Log2 RPKM of miRNA") + scale_y_discrete()
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

make_cor_plots_median_log2 <- function(chi, title){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=mad(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing_Samples=is_complete(value))
  ggplot(tmp) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing_Samples)) +
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Median(RPKM)", y="MAD(RPKM)") + xlim(c(0,20)) + ggtitle(title)#+ geom_abline(aes(color="red"))
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
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Median(RPKM)", y="MAD(RPKM)") + geom_abline()
  #     scale_size_continuous(name="MAD(RPKM)", breaks=c(25000, 50000, 75000 ), labels=c("25,000", "50,000", "75,000"))
}

isComRPKM <- function(chi){
  return (!(0.0 %in% chi))
}


make_cor_plots_mean <- function(chi, title){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=mean(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=sd(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing=isComRPKM(value))
  ggplot(tmp) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing)) +
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Mean(RPKM)", y="SD(RPKM)") + geom_abline() + ggtitle(title)
}

make_cor_plots_mean_zoom <- function(chi){
  tmp <- chi %>% gather(gene, "value")
  colnames(tmp)[2] <- "sample"
  tmp <- tmp %>% group_by(gene) %>% mutate(avg=mean(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(deviation=sd(value, na.rm=T))
  tmp <- tmp %>% group_by(gene) %>% mutate(No_Missing=isComRPKM(value))
  ggplot(tmp %>% group_by(gene)) + geom_point(aes(x=avg, y=deviation, size=deviation, color=No_Missing)) +
    theme(axis.title.y = element_text(vjust=-0.001)) +
    labs(x="Mean(RPKM)", y="SD(RPKM)") + geom_abline() + xlim(0, 100) + ylim(c(0, 100))
}


ind_kde_plot <- function(x, k=.10){
  y <- x %>% sample_frac(k) %>% gather(gene)
  ret <- ggplot(y, aes(x = value, fill = gene)) + geom_density(alpha = .2) + labs(x = "Log2 RPKM") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0)
  return (ret)
}


make_KDE_plots <- function(chi, title, k=100, individual=F){
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
    
    ggplot(y, aes(vals)) + geom_density(alpha = 0.4)  + ggtitle(title) + labs(x = "Log(2) RPKM") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0)   
  }
  else{
    z_tot <<- unlist(zed)
#     print(ggplot(zed, aes(x = zed)) + geom_density(alpha = 0.4) + scale_fill_hue(l=30) + ggtitle(title) + labs(x = "RPKM (Log2") + guides(fill = F) + xlim(-5, 20) + ylim(0.0, 1.0))

  }
}


## log2 Expression, for tri plots
ACC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ACC.RPKM.log2.munged.txt", sep="\t"))
ACC.RPKM.log2.munged <- remove_non_tumor(ACC.RPKM.log2.munged)
make_tri_plot(ACC.RPKM.log2.munged, "ACC", c(0, 150))

BLCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BLCA.RPKM.log2.munged.txt", sep="\t"))
BLCA.RPKM.log2.munged <- remove_non_tumor(BLCA.RPKM.log2.munged)
make_tri_plot(BLCA.RPKM.log2.munged, "BLCA", c(0, 400))

BRCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.log2.munged.txt", sep="\t"))
BRCA.RPKM.log2.munged <- remove_non_tumor(BRCA.RPKM.log2.munged)
make_tri_plot(BRCA.RPKM.log2.munged, "BRCA", c(0, 800))

CESC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/CESC.RPKM.log2.munged.txt", sep="\t"))
CESC.RPKM.log2.munged <- remove_non_tumor(CESC.RPKM.log2.munged)
make_tri_plot(CESC.RPKM.log2.munged, "CESC", c(0, 300))

COADREAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COADREAD.RPKM.log2.munged.txt", sep="\t"))
COADREAD.RPKM.log2.munged <- remove_non_tumor(COADREAD.RPKM.log2.munged)
make_tri_plot(COADREAD.RPKM.log2.munged, "COADREAD", c(0, 400))

COAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COAD.RPKM.log2.munged.txt", sep="\t"))
COAD.RPKM.log2.munged <- remove_non_tumor(COAD.RPKM.log2.munged)
make_tri_plot(COAD.RPKM.log2.munged, "COAD", c(0, 400))

DLBC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/DLBC.RPKM.log2.munged.txt", sep="\t"))
DLBC.RPKM.log2.munged <- remove_non_tumor(DLBC.RPKM.log2.munged)
make_tri_plot(DLBC.RPKM.log2.munged, "DLBC", c(0, 100))

ESCA.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ESCA.RPKM.log2.munged.txt", sep="\t"))
ESCA.RPKM.log2.munged <- remove_non_tumor(ESCA.RPKM.log2.munged)
make_tri_plot(ESCA.RPKM.log2.munged, "ESCA", c(0, 200))

HNSC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/HNSC.RPKM.log2.munged.txt", sep="\t"))
HNSC.RPKM.log2.munged <- remove_non_tumor(HNSC.RPKM.log2.munged)
make_tri_plot(HNSC.RPKM.log2.munged, "HNSC", c(0, 600))

KICH.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KICH.RPKM.log2.munged.txt", sep="\t"))
KICH.RPKM.log2.munged <- remove_non_tumor(KICH.RPKM.log2.munged)
make_tri_plot(KICH.RPKM.log2.munged, "KICH", c(0, 200))

KIRC.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KIRC.RPKM.log2.munged.txt", sep="\t")) 
KIRC.RPKM.log2.munged <- remove_non_tumor(KIRC.RPKM.log2.munged)
make_tri_plot(KIRC.RPKM.log2.munged, "KIRC", c(0, 400))

LAML.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/LAML.RPKM.log2.munged.txt", sep="\t"))
make_tri_plot(LAML.RPKM.log2.munged, "LAML", c(0, 200))

PAAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PAAD.RPKM.log2.munged.txt", sep="\t"))
PAAD.RPKM.log2.munged <- remove_non_tumor(PAAD.RPKM.log2.munged)
make_tri_plot(PAAD.RPKM.log2.munged, "PAAD", c(0, 200))

PRAD.RPKM.log2.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PRAD.RPKM.log2.munged.txt", sep="\t"))
PRAD.RPKM.log2.munged <- remove_non_tumor(PRAD.RPKM.log2.munged)
make_tri_plot(PRAD.RPKM.log2.munged, "PRAD", c(0, 400))


all_log2_frames <- list(
  ACC.RPKM.log2.munged, BLCA.RPKM.log2.munged, BRCA.RPKM.log2.munged,
COAD.RPKM.log2.munged, COADREAD.RPKM.log2.munged, CESC.RPKM.log2.munged, DLBC.RPKM.log2.munged,
ESCA.RPKM.log2.munged, HNSC.RPKM.log2.munged, KICH.RPKM.log2.munged, KIRC.RPKM.log2.munged, LAML.RPKM.log2.munged,
PAAD.RPKM.log2.munged, PRAD.RPKM.log2.munged)

all_frames <- list(
  ACC.RPKM.munged, BLCA.RPKM.munged, BRCA.RPKM.munged,
  COAD.RPKM.munged, COADREAD.RPKM.munged, CESC.RPKM.munged, DLBC.RPKM.munged,
  ESCA.RPKM.munged, HNSC.RPKM.munged, KICH.RPKM.munged, KIRC.RPKM.munged, LAML.RPKM.munged,
  PAAD.RPKM.munged, PRAD.RPKM.munged)

log_names <- list("ACC", "BLCA", "BRCA", "COAD", "COADREAD", "CESC", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD")

for (i in c(1:14)){
#   cat(str(i))
#   cat("\n")
#   colnames(all_log2_frames[i])[1] <- "gene"
  message(log_names[i])
  print(make_cor_plots_mean_log2(all_log2_frames[[i]], log_names[[i]]))
}

for (i in c(1:14)){
  #   cat(str(i))
  #   cat("\n")
  #   colnames(all_log2_frames[i])[1] <- "gene"
  message(log_names[i])
  print(make_L_plots_dplyr(all_log2_frames[[i]], log_names[[i]]))
}

for (i in c(1:14)){
  #   cat(str(i))
  #   cat("\n")
  colnames(all_frames[[i]])[1] <- "gene"
  message(log_names[i])
  print(make_cor_plots_mean(all_frames[[i]], log_names[[i]]))
}

##
##Overall KDEs
##

## KDE of Total Dists
brca <- as_data_frame(data.frame(vals=unlist(BRCA.RPKM.log2.munged), sample="BRCA"))
ggplot(brca, aes(x=vals)) + geom_density(alpha=.4) + xlim(-5, 20) + ggtitle("BRCA Overall Distribution") + xlab("Log2 RPKM") + theme(axis.title.y = element_text(vjust=-0.001))


## split (HE/LE KDE)
is_complete <- function(x){
  return (countNAs(x) == 0)
}

is_high <- function(x, lim=5.0){
  return (!(mean(x) <= lim | is.na(mean(x, na.rm=T))))
}
acc <- gather(ACC.RPKM.log2.munged, gene, "value") 
colnames(acc)[2] <- "sample"
acc <- acc %>% group_by(gene) %>% mutate(isCom=is_complete(value))
acc <- acc %>% group_by(gene) %>% mutate(isHigh=is_high(value, 4.0))
ggplot(acc, aes(x=value, fill=isCom)) + geom_density(alpha=.3)
ggplot(acc, aes(x=value, fill=isHigh)) + geom_density(alpha=.3)

brca <- gather(BRCA.RPKM.log2.munged, gene, "value")
colnames(brca)[2] <- "sample"
brca <- brca %>% group_by(gene) %>% mutate(isCom=is_complete(value))
brca <- brca %>% group_by(gene) %>% mutate(isHigh=is_high(value, 4.0))
ggplot(brca, aes(x=value, fill=isHigh)) + geom_density(alpha=.4)  + xlim(c(-5, 20)) + ylim(c(0,.5)) + scale_fill_discrete(name="AVG Log2(RPKM)", breaks=c("TRUE", "FALSE"), labels=c(">5 RPKM", "<=5 RPKM")) + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for BRCA")
ggplot(brca, aes(x=value, fill=isCom)) + geom_density(alpha=.4)  + xlim(c(-5, 20)) + ylim(c(0,.5)) + scale_fill_discrete(name="Missing Values", breaks=c("TRUE", "FALSE"), labels=c("None", "At least one sample missing")) + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for BRCA")



ggplot(brca) +
  geom_density(alpha=.4, aes(x=value, fill=isHigh))+
  xlim(c(-5, 20)) + ylim(0,.25) +
  geom_density(alpha=.04, aes(x=value, color=isCom)) +
  ggtitle("Comparing Mean RPKM\n and Missing Values") +
  theme(axis.title.y = element_text(vjust=-0.001))
#scale_colour_discrete(name="Missing Values", breaks=c("TRUE", "FALSE"), labels=c("None", "At least one sample missing")) +
#scale_fill_discrete(name="Mean RPKM of miRNA", breaks=c("True", "False"), labels=c(">= 4 RPKM", "<4 RPKM"))

acc$c_type <- "ACC"
brca$c_type <- "BRCA"
combi <- bind_rows(acc, brca)
ggplot(combi, aes(x=value, fill=c_type)) + geom_density(alpha=.4) + scale_fill_discrete(name="Cancer Type") + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for ACC and BRCA")


brca_full <- as_data_frame(data.frame(vals=unlist(BRCA.RPKM.log2.munged), sample="BRCA"))

ggplot(brca, aes(x=value)) + 
  geom_density(aes(x=value, colour="green"), alpha=0.4) +
  geom_density(aes(x=value, fill=isCom), alpha=0.2) +
  xlim(-5,20) + ggtitle("KDE for BRCA") + ylim(0, 1.0)+
theme(axis.title.y = element_text(vjust=-0.001)) + theme(legend.position="none")
legend("topright", c("Overall", "miRNAs w/ Missing Samples", "miRNAs In All Samples"), col=c(2:4), pch=19, cex=0.4)

ggplot(brca, aes(x=value, fill=gene)) + geom_density(alpha=.4)  + xlim(c(-5, 20)) + ylim(c(0,1.0))  + theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("KDE for BRCA, miRNA") + theme(legend.position="none")


## Look for microRNAs with missing values and high means



is_rpkm_complete <- function(x){
  return (!(0.0 %in% x))
}

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

##
## Raw Expression dendrograms
## for both LE and HE expressed genes.
ACC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ACC.RPKM.munged.txt", sep="\t"))
acc <- prep_for_dendro(ACC.RPKM.munged, 1)
  
BLCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BLCA.RPKM.munged.txt", sep="\t"))
blca <- prep_for_dendro(BLCA.RPKM.munged, 2)

BRCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.munged.txt", sep="\t"))
brca <- prep_for_dendro(BRCA.RPKM.munged, 3)

CESC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/CESC.RPKM.munged.txt", sep="\t"))
cesc <- prep_for_dendro(CESC.RPKM.munged, 4)

COADREAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COADREAD.RPKM.munged.txt", sep="\t"))
coadread <- prep_for_dendro(COADREAD.RPKM.munged, 5)

COAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/COAD.RPKM.munged.txt", sep="\t"))
coad <- prep_for_dendro(COAD.RPKM.munged, 6)

DLBC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/DLBC.RPKM.munged.txt", sep="\t"))
dlbc <- prep_for_dendro(DLBC.RPKM.munged, 7)

ESCA.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/ESCA.RPKM.munged.txt", sep="\t"))
esca <- prep_for_dendro(ESCA.RPKM.munged, 8)

HNSC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/HNSC.RPKM.munged.txt", sep="\t"))
hnsc <- prep_for_dendro(HNSC.RPKM.munged, 9)

KICH.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KICH.RPKM.munged.txt", sep="\t"))
kich <- prep_for_dendro(KICH.RPKM.munged, 10)

KIRC.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/KIRC.RPKM.munged.txt", sep="\t"))
kirc <- prep_for_dendro(KIRC.RPKM.munged, 11)

LAML.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/LAML.RPKM.munged.txt", sep="\t"))
laml <- prep_for_dendro(LAML.RPKM.munged, 12)

PAAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PAAD.RPKM.munged.txt", sep="\t"))
paad <- prep_for_dendro(PAAD.RPKM.munged, 13)

PRAD.RPKM.munged <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/PRAD.RPKM.munged.txt", sep="\t"))
prad <- prep_for_dendro(PRAD.RPKM.munged, 14)

tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))
clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab="")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.35)

##
## Dendros for just LE genes
##
prep_for_dendroLE<- function(chi, type_val, le=T){
  chi <- chi %>% gather(sample)
  colnames(chi)[2] <- "gene"
  chi <- chi %>% group_by(gene) %>% mutate(complete=is_rpkm_complete(value))
  if (le == T){
    chi <- chi %>% filter(complete == F)
  }
  else{
    chi <- chi %>% filter(complete == T)
  }
  chi <- chi %>% select(-complete) %>% spread(gene, value)
  chi <- chi %>% mutate(type=type_val)
  return (chi)
}
acc <- prep_for_dendroLE(acc, 1)
blca <- prep_for_dendroLE(blca, 2)
brca <- prep_for_dendroLE(brca, 3)
cesc <- prep_for_dendroLE(cesc, 4)
coadread <- prep_for_dendroLE(coadread, 5)
coad <- prep_for_dendroLE(coad, 6)
dlbc <- prep_for_dendroLE(dlbc, 7)
esca <- prep_for_dendroLE(esca, 8)
hnsc <- prep_for_dendroLE(hnsc, 9)
kich <- prep_for_dendroLE(kich, 10)
kirc <- prep_for_dendroLE(kirc, 11)
laml <- prep_for_dendroLE(laml, 12)
paad <- prep_for_dendroLE(paad, 13)
prad <- prep_for_dendroLE(prad, 14)


tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))
clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab="")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.35)


##
## Dendros for HE
##
acc <- prep_for_dendroLE(acc, 1, F)
blca <- prep_for_dendroLE(blca, 2, F)
brca <- prep_for_dendroLE(brca, 3, F)
cesc <- prep_for_dendroLE(cesc, 4, F)
coadread <- prep_for_dendroLE(coadread, 5, F)
coad <- prep_for_dendroLE(coad, 6, F)
dlbc <- prep_for_dendroLE(dlbc, 7, F)
esca <- prep_for_dendroLE(esca, 8, F)
hnsc <- prep_for_dendroLE(hnsc, 9, F)
kich <- prep_for_dendroLE(kich, 10, F)
kirc <- prep_for_dendroLE(kirc, 11, F)
laml <- prep_for_dendroLE(laml, 12, F)
paad <- prep_for_dendroLE(paad, 13, F)
prad <- prep_for_dendroLE(prad, 14, F)

tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))
tmp_fin[is.na(tmp_fin)] <- 0.0
clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("sample"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000, xlab="")
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.35)


##
## A plot of number of missing values on mean RPKM
## shows that above 4-6RPKM, there are no more missing values.
##

make_inverted_L_plot <- function(chi){
#   x <- x %>% select(contains("TCGA"))
  chi <- chi %>% gather(gene)
  colnames(x)[2] <- "sample"
  chi <- chi %>% group_by(gene) %>% mutate(avg=mean(value, na.rm=T))
  chi <- chi %>% group_by(gene) %>% mutate(numNAs=countNAs(value))
  ggplot(chi) + geom_point(aes(x=mean, y=noNAs))
}
BRCA.RPKM.log2.munged.transposed <- as_data_frame(read.csv("/home/eric/sync/data/miRNA-seq/BRCA.RPKM.log2.munged.txt.transposed", sep="\t"))
BRCA.RPKM.log2.munged.transposed_no_gene <- BRCA.RPKM.log2.munged.transposed %>% select(-contains("gene"))

x <- apply(BRCA.RPKM.log2.munged.transposed_no_gene, 2, mean, na.rm=T)
n <- colnames(BRCA.RPKM.log2.munged.transposed_no_gene)
y <- apply(BRCA.RPKM.log2.munged.transposed_no_gene, 2, countNAs)

mac <- as_data_frame(data.frame(nam = n, means = x, numNAs = y))
mac %>% arrange(desc(means)) -> mac
qplot(mac$means, mac$numNAs, xlim = c(-5, 20)) +
  xlab("Mean Log2 RPKM") + ylab("Number of Samples Missing miRNA") +
  theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("Number of Missing Values Decreases with Increasing Mean RPKM (BRCA)")

qplot(mac$means, mac$numNAs, xlim = c(-5, 20), ylim=c(0,5)) +
  xlab("Mean Log2 RPKM") + ylab("Number of Samples Missing miRNA") +
  theme(axis.title.y = element_text(vjust=-0.001)) + ggtitle("(BRCA)")






##
## Dendrograms, but with only LE genes
acc <- ACC.RPKM.log2.munged %>% gather(gene, "value") 
colnames(acc)[2] <- "sample"
acc <- acc %>% group_by(gene) %>% mutate(complete=is_complete(value)) %>% filter(complete==T)
acc <- acc %>% spread(gene, value) %>% select(-complete) %>%

tmp_fin <- do.call(bind_rows, list(acc, blca, brca,
                                   cesc, coad, coadread, dlbc, esca,
                                   hnsc, kich, kirc, laml, paad, prad))

clus <- tmp_fin %>% select(-contains("type")) %>% select(-contains("HYBRID"))
hh <- hclust(dist(clus))
ColorDendrogram(hh, y=as.numeric(tmp_fin$type), branchlength = 300000)
legend("topright", c("ACC", "BLCA", "BRCA", "CESC", "COAD", "COADREAD", "DLBC", "ESCA", "HNSC", "KICH", "KIRC", "LAML", "PAAD", "PRAD"), col=c(2:15), pch=19, cex=0.35)


log_mad <- function(x, na.rm = F){
  ret <- log(mad(x, na.rm))
  return(ret)
}
x <- apply(ACC.RPKM.log2.munged.transposed[,-1], 2, mean, na.rm=T)
y <- apply(ACC.RPKM.log2.munged.transposed[,-1], 2, log_mad, na.rm=T)
qplot(x, y, main = "Log(MAD) on Log(Mean) of Expression")



## Normal cell analysis
norm_breast <- as_data_frame(read.csv("/home/eric/htseq_out_mirna_ERR86_RPKM.txt", sep="\t"))

## BRCA dist fitting
# fit_distribs_helper <- function(d){
#   ## A fitting function borrowed from World of Piggy
#   tmp <- as.numeric(d)
#   tmp <- tmp[!is.na(tmp)]
#   ret <- c()
#   for (i in c("normal", "gamma")){
#     fit <- fitdistr(tmp, i)
#     est_shape <- fit$estimate[[1]]
#     est_rate <- fit$estimate[[2]]
#     ks = ks.test(tmp, "pgamma", shape=est_shape, rate=est_rate)
#     
#   }
# }


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
detach("package:dplyr", unload=T)
require(MASS)


BRCA.mRNA <- as_data_frame(read.csv("/home/eric/sync/data/brca/BRCA.mRNAseq_RPKM_log2.txt", sep="\t"))
colnames(BRCA.mRNA)[1] <- "gene"
BRCA.mRNA <- remove_non_tumor(BRCA.mRNA)
BRCA.mRNA <- BRCA.mRNA %>% gather(gene)
colnames(BRCA.mRNA)[2] <- "sample"
BRCA.mRNA = BRCA.mRNA %>% group_by(gene) %>% mutate(isComplete=is_complete(value))

ggplot(BRCA.mRNA %>% sample_frac(.05)) + geom_density(aes(x=value, position="fill")) +
 geom_density(aes(x=value, color=isComplete, position="fill"))

qplot(value, ..count.., data=BRCA.mRNA %>% sample_frac(.50), geom="density", color=isComplete, adjust=1, xlim=c(-5, 20))


brca <- BRCA.RPKM.log2.munged %>% gather(gene, value)
colnames(brca)[2] <- "sample"
brca <- brca %>% group_by(gene) %>% mutate(is_com=is_complete(value))
brca <- brca %>% group_by(gene) %>% mutate(isH=is_high(value))
brca <- brca %>% group_by(gene) %>% mutate(med=median(value, na.rm=T))
brca <- brca %>% group_by(gene) %>% mutate(mm=mad(value, na.rm=T))
brca <- brca %>% group_by(gene) %>% mutate(numNAs=countNAs(value))


qplot(value, ..count.., data=brca %>% sample_frac(1.0), geom="density", color=isComplete, adjust=1, xlim=c(-5, 20))

get_lost_expr <- function(chi){
  chi %>% gather(gene, "value") -> chi
  colnames(chi)[2] <- "sample"
  chi %>% group_by(gene) %>% mutate(avg=median(value, na.rm=T)) -> chi
  chi %>% group_by(gene) %>% mutate(numNAs=countNAs(value)) -> chi
  ret <- chi %>% group_by(gene) %>% filter(avg > 5.0) %>% filter(numNAs > 0)
  return (ret)
}


for (i in c(1:14)){
  #   cat(str(i))
  #   cat("\n")
  #   colnames(all_log2_frames[i])[1] <- "gene"
#   message(log_names[i])
  tmp <- get_lost_expr(all_log2_frames[[i]])
  ret <- spread(tmp, sample, value)
  message(log_names[i])
  message(paste(ret$gene, sep=" ", collapse=" "))
}

## Figure 1

## Figure 2

## Figure 3

## Figure 4

## Figure 5

## Figure 6

## Figure 7

## Figure 8

## Figure 9

## Figure 10