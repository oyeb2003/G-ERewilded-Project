library(ggplot2)
library(circlize)
library(gridExtra)
library(reshape)
library(ape)
library(scales)
library(ggrepel)
library(MDMR)
library(ggsci)
library(ggsignif)

#This code and/or a modification of it was used for ploting all the PCA plots in this manuscript

setwd("/Volumes/oyesolaoo$/My Documents/Rewilding Project/Flow Cytometry/Blood Analysis Unsupervised")

metadata = read.csv("FACS_Blood_Unsupervised_Labelled1.csv",header = T, row.names = 1)
head(metadata)
colnames(metadata)

input = metadata[,c(7:21)]
input 
metadata[, 2:4] <- lapply(metadata[,2:4], as.factor)
effectors= metadata[,2:4]
effectors

  D <- dist(input, method = "euclidean")
  res <- pcoa(D)
  plotter <- data.frame(res$vectors, col_var=metadata$Genotype, shape_var=metadata$Environment)
  
  
  
  #good_plotter <- subset(plotter, Axis.1 > quantile(plotter$Axis.1, 0.05) &
  #                         Axis.1 < quantile(plotter$Axis.1, 0.95))
  #good_plotter <- subset(good_plotter, Axis.2 > quantile(plotter$Axis.2, 0.05) &
  #                         Axis.2 < quantile(plotter$Axis.2, 0.95))
  
  
  g0=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape = shape_var)) + 
    geom_point(size=8) +
    guides(color=guide_legend(title="Genotype"), 
           shape=guide_legend(title="Environment")) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_shape_manual(values = c(18, 3)) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) + 
    #scale_size_manual(values = c(2,3,10))+
    ggtitle("Genotype")+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black",face='bold',hjust = 0.5))
  g0
  
  pr.coo=res$vectors
  plot.axes=c(1,2)
  n <- nrow(input)
  points.stand <- scale(pr.coo[,plot.axes])
  S <- cov(input, points.stand)
  U <- S %*% diag((res$values$Eigenvalues[plot.axes]/(n-1))^(-0.5))
  colnames(U) <- colnames(pr.coo[,plot.axes])
  
  PC = res
  data <- data.frame(obsnames=row.names(PC$vectors), PC$vectors[,1:2])
  datapc <- data.frame(varnames=rownames(U), U*100)
  datapc$var1 <- rescale(datapc$Axis.1, c(min(data$Axis.1),max(data$Axis.1)))
  datapc$var2 <- rescale(datapc$Axis.1, c(min(data$Axis.2),max(data$Axis.2)))
  
  datapc$mult <- abs(datapc$Axis.1*datapc$Axis.2)
  datapc <- datapc[order(datapc$mult, decreasing = T),]
  datapc2 = datapc
  datapc2 = datapc[1:12,]
  
  g_seg1=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var)) + 
    geom_point(size=5) +
    guides(color=guide_legend(title="Environment"))+
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    ggtitle("Environment")+
    theme_bw() + theme(legend.position = "none") +
    theme(plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5)) 
  g_seg1
  
  g_seg2=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape = shape_var)) + 
    geom_point(size=8) +
    guides(color=guide_legend(title="Genotype"), 
           shape=guide_legend(title="Environment")) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_shape_manual(values = c(18, 3)) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) + 
    scale_size_manual(values = c(2,3,6))+
    ggtitle("Genotype")+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5))
  g_seg2
  
  g_seg=ggplot(plotter, aes(Axis.1, Axis.2)) +
    geom_point(size=20, alpha=0) + ggtitle("Loadings for Clusters") +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    geom_segment(data=datapc2, aes(x=0, y=0, xend=Axis.1, yend=Axis.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
    geom_label_repel(data=datapc2, aes(x=Axis.1, y=Axis.2, label=varnames), 
                     size = 6, force=10, segment.alpha=0.5) +
  theme_void()
  g_seg
  
  
  pc1_hist <- ggplot(plotter, aes(Axis.1, color = shape_var, fill= shape_var)) +
    geom_density(alpha=0.6) + 
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    scale_fill_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    theme_void() + theme(legend.position = "none")
  pc1_hist
    
  pc2_hist <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    coord_flip() +
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3", "forestgreen")) +
    scale_fill_manual(values = c("mediumpurple1", "red3", "forestgreen")) +
    theme_void() + theme(legend.position = "none") 
  pc2_hist
  
  #pc_fake <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    #geom_density(alpha=0.6) + 
    #coord_flip() +
    #scale_y_reverse() + 
    #scale_color_manual(values = c("white", "white", "white")) +
    #scale_fill_manual(values = c("white", "white", "white")) +
    #theme_void() + theme(legend.position = "none") 
  #pcfake


  PCA Analysis and Code MLN Flow Cytometry data 

  library(ggplot2)
library(circlize)
library(gridExtra)
library(reshape)
library(ape)
library(scales)
library(ggrepel)
library(MDMR)
library(ggsci)
library(ggsignif)


setwd("/Volumes/oyesolaoo$/My Documents/Rewilding Project/Flow Cytometry/MLN Analysis Unsupervised")

metadata = read.csv("MLN_ICS_Compositional data_named.csv",header = T, row.names = 1)
head(metadata)

input = metadata[,c(4:18)]
input
metadata[, 1:3] <- lapply(metadata[,1:3], as.factor)
effectors= metadata[,1:3]
effectors

  D <- dist(input, method = "euclidean")
  res <- pcoa(D)
  print(res)
  class(res)
  str(res)
  
  sum <- as.data.frame(res$vectors)
  sum
  write.csv(sum, "PCA_April.csv")
  
  plotter <- data.frame(res$vectors, col_var=metadata$Strain, shape_var=metadata$Location)
  
  
  
  #good_plotter <- subset(plotter, Axis.1 > quantile(plotter$Axis.1, 0.05) &
  #                         Axis.1 < quantile(plotter$Axis.1, 0.95))
  #good_plotter <- subset(good_plotter, Axis.2 > quantile(plotter$Axis.2, 0.05) &
  #                         Axis.2 < quantile(plotter$Axis.2, 0.95))
  
  
  g0=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape = shape_var)) + 
    geom_point(size=8) +
    guides(color=guide_legend(title="Genotype"), 
           shape=guide_legend(title="Environment")) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_shape_manual(values = c(18, 3)) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) + 
    #scale_size_manual(values = c(2,3,10))+
    ggtitle("Genotype")+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black",face='bold',hjust = 0.5))
  g0
  
  pr.coo=res$vectors
  plot.axes=c(1,2)
  n <- nrow(input)
  points.stand <- scale(pr.coo[,plot.axes])
  S <- cov(input, points.stand)
  U <- S %*% diag((res$values$Eigenvalues[plot.axes]/(n-1))^(-0.5))
  colnames(U) <- colnames(pr.coo[,plot.axes])
  
  PC = res
  data <- data.frame(obsnames=row.names(PC$vectors), PC$vectors[,1:2])
  datapc <- data.frame(varnames=rownames(U), U*100)
  datapc$var1 <- rescale(datapc$Axis.1, c(min(data$Axis.1),max(data$Axis.1)))
  datapc$var2 <- rescale(datapc$Axis.1, c(min(data$Axis.2),max(data$Axis.2)))
  
  datapc$mult <- abs(datapc$Axis.1*datapc$Axis.2)
  datapc <- datapc[order(datapc$mult, decreasing = T),]
  datapc2 = datapc
  datapc2 = datapc[1:12,]
  
  g_seg1=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var)) + 
    geom_point(size=5) +
    guides(color=guide_legend(title="Environment"))+
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    ggtitle("Environment")+
    theme_bw() + theme(legend.position = "none") +
    theme(plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5)) 
  g_seg1
  
  g_seg2=ggplot(plotter, aes(Axis.1, Axis.2, color=col_var, shape = shape_var)) + 
    geom_point(size=8) +
    guides(color=guide_legend(title="Genotype"), 
           shape=guide_legend(title="Environment")) +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    scale_shape_manual(values = c(18, 3)) +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) + 
    scale_size_manual(values = c(2,3,6))+
    ggtitle("Genotype")+
    theme_bw() +
    theme(axis.text = element_text(size=12, color="black"),
          axis.title = element_text(size=13, color="black"),
          legend.text = element_text(size=12, color="black"),
          legend.title = element_text(size=13, color="black"),
          plot.title = element_text(size=14, color="black", 
                                    face='bold',hjust = 0.5))
  g_seg2
  
  g_seg=ggplot(plotter, aes(Axis.1, Axis.2)) +
    geom_point(size=3, alpha=0) + ggtitle("Loadings for Clusters") +
    ylab(paste0("PCo2 ", round(res$values$Broken_stick[2]*100,2), "% expl. variation")) +
    xlab(paste0("PCo1 ", round(res$values$Broken_stick[1]*100,2), "% expl. variation")) +
    geom_segment(data=datapc2, aes(x=0, y=0, xend=Axis.1, yend=Axis.2), 
                 arrow=arrow(length=unit(0.2,"cm")), alpha=0.5) +
    geom_label_repel(data=datapc2, aes(x=Axis.1, y=Axis.2, label=varnames), 
                     size = 4, force=4, segment.alpha=0.5) +
    theme_void() 
  g_seg
  
  
  pc1_hist <- ggplot(plotter, aes(Axis.1, color =  col_var, fill= col_var)) +
    geom_density(alpha=0.6) + 
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    scale_fill_manual(values = c("mediumpurple1", "red3","forestgreen")) +
    theme_void() + theme(legend.position = "none")
  pc1_hist
    
  pc2_hist <- ggplot(plotter, aes(Axis.2, color = shape_var, fill= shape_var)) +
    geom_density(alpha=0.6) + 
    coord_flip() +
    scale_y_reverse() +
    scale_color_manual(values = c("mediumpurple1", "red3", "forestgreen")) +
    scale_fill_manual(values = c("mediumpurple1", "red3", "forestgreen")) +
    theme_void() + theme(legend.position = "none") 
  pc2_hist
  
  #pc_fake <- ggplot(plotter, aes(Axis.2, color = col_var, fill= col_var)) +
    #geom_density(alpha=0.6) + 
    #coord_flip() +
    #scale_y_reverse() + 
    #scale_color_manual(values = c("white", "white", "white")) +
    #scale_fill_manual(values = c("white", "white", "white")) +
    #theme_void() + theme(legend.position = "none") 
  #pcfake
  
combine <- grid.arrange(pc2_hist, g_seg1, pc_fake, pc1_hist, ncol=2)

ggsave("Combine_plot_environment.pdf", combine)


  
combine <- grid.arrange(pc2_hist, g_seg1, pc_fake, pc1_hist, ncol=2)

ggsave("Combine_plot_environment.pdf", combine)

