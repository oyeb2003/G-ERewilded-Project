# G-ERewilded-Project
All codes used for the Genetic by environment interaction analysis and differential gene expression analysis performed in the manuscript: Genetic and Environmental interactions contribute to immune variation in rewilded mice by Oyesola et al. 2023

If  you need more information about this - please contact - Dr. Andrea Graham

setwd("~/Documents/live DATA/Stony Ford/2021 for R")
library(tidyverse)
library(MASS)
Fig_8C<-read.csv("Worms_scRNA.csv",sep=",",dec=".",header=T)
Fig_8C$Genotype<-as.factor(Fig_5C$Genotype)
Fig_8C$Genotype<-relevel(Fig_5C$Genotype,ref = "C57")
Worms_RNA_PC2<-glm.nb(Worm.count~Genotype+Location+PC2,data=Fig_8C)
summary(Worms_RNA_PC2)
drop1(Worms_RNA_PC2,test="Chisq")

model_sumPC2 <- summary(Worms_RNA_PC2)
estimatesPC2_1 <- rnorm(1000,
                       mean=model_sumPC2$coefficients["(Intercept)","Estimate"],
                       sd=model_sumPC2$coefficients["(Intercept)","Std. Error"])
estimatesPC2_2 <- rnorm(1000,
                       mean=model_sumPC2$coefficients["Genotype129","Estimate"],
                       sd=model_sumPC2$coefficients["Genotype129","Std. Error"])
estimatesPC2_3 <- rnorm(1000,
                       mean=model_sumPC2$coefficients["GenotypePWK","Estimate"],
                       sd=model_sumPC2$coefficients["GenotypePWK","Std. Error"])
estimatesPC2_4 <- rnorm(1000,
                       mean=model_sumPC2$coefficients["LocationSF","Estimate"],
                       sd=model_sumPC2$coefficients["LocationSF","Std. Error"])
estimatesPC2_5 <- rnorm(1000,
                       mean=model_sumPC2$coefficients["PC2","Estimate"],
                       sd=model_sumPC2$coefficients["PC2","Std. Error"])
estimatesPC2 <- as.data.frame(cbind(estimatesPC2_1,estimatesPC2_2,estimatesPC2_3,estimatesPC2_4,estimatesPC2_5))
estimatesPC2 <- tidyr::pivot_longer(estimatesPC2,1:5,names_to="Parameter")

ggplot(estimatesPC2,aes(value,Parameter))+
  geom_violin(aes(fill=Parameter),scale="width")+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text.y=element_text(angle=45,hjust=0.5,size=10),
        axis.line=element_line(lineend="round"))+
  labs(y="Parameter in model",
       x="Model-estimated effect on worm burden")+
  scale_y_discrete(labels=c("C57 indoors (intercept)",
                            "129 vs. C57\nstrain","PWK vs. C57\nstrain",
                            "Outdoors","PC2"))+
  scale_fill_manual(values=MetBrewer::met.brewer("Hokusai2",5,type="discrete"))+
  guides(fill="none")




