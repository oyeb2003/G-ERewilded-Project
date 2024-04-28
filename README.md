# G-ERewilded-Project
##Effect size calculation of Spectral Cytometry data in Figure 1 and Figure 2

#Installing and loading packages
library(MDMR)
#install.packages("lme4")
library(lme4)
library(tidyverse)
#install.packages("broom")
library(broom)

#choosing the explanatory variables 
setwd("/Volumes/oyesolaoo$/My Documents/Rewilding Project/Block 2/Flow Data/Blood/2022/Analysis MDMR")
X.mdmr = read.csv("x_bt.mdmr2.csv",header = T, row.names = 1)
X.mdmr
#X.mdmr$Block <- as.factor(X.mdmr$Block)
#x.mdmr2 = read.csv("Block1_Block_2_Mouse_Metadata_ELISA_3_24.csv",header = T, row.names = 1)


#returns the dimension (e.g. the number of columns and rows) of a matrix, array or data frame
dim(X.mdmr)

#choosing the response variable 
Y.mdmr = read.csv("y_bt.mdmr.csv",header = T, row.names = 1)

#returns the dimension (e.g. the number of columns and rows) of a matrix, array or data frame
dim(Y.mdmr)

#All variables are standardized to mean zero and unit variance.

Y.mdmr <- scale(Y.mdmr,center = TRUE, scale = TRUE)

# Compute distance matrix for the y variable (quantitative/response variable)
D <- dist(Y.mdmr, method = "euclidean")

#A Conduct MDMR (regress the quantitative/response variable unto the explanatory variables)
#Option A: Passing the Distance Matrix - note the perm.p = False function
mdmr.res <- mdmr(X = X.mdmr, D = D)


# Check results
#summary(mdmr.res) 

sum <- as.data.frame(summary(mdmr.res))

write.csv(sum, "summary_results.csv")

#Accounting for grouping variable

tempframe <- model.frame(~ Env*Infection*Geno, data = X.mdmr)
temp.mdmr <- model.matrix(~ Env*Infection*Geno, data = X.mdmr)
tempframe
temp.mdmr

res <- mdmr(X = temp.mdmr, D = D)

sum_interaction <- as.data.frame(summary(res))
sum_interaction
write.csv(sum_interaction, "summary_interaction.csv")
![image](https://github.com/oyeb2003/G-ERewilded-Project/assets/79373782/f022b07b-6a33-403f-a980-91d112d722f1)

