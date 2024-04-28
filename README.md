# G-ERewilded-Project
All codes used for the Genetic by environment interaction analysis and differential gene expression analysis performed in the manuscript: Genetic and Environmental interactions contribute to immune variation in rewilded mice by Oyesola et al. 2023
#Installing and loading packages
library(MDMR)
#install.packages("lme4")
library(lme4)

setwd("/Volumes/oyesolaoo$/My Documents/Rewilding Project/ELISA/Serum/Analysis/Block2")
#choosing the explanatory variables 
X.mdmr = read.csv("x_bt.mdmr.csv",header = T, row.names = 1)
X.mdmr
table(X.mdmr$Genotype)
table(X.mdmr$Location)
table(X.mdmr$Infection_Status)
table(X.mdmr$Block)
table(X.mdmr$Genotype, X.mdmr$Location, X.mdmr$Infection_Status)



#returns the dimension (e.g. the number of columns and rows) of a matrix, array or data frame
dim(X.mdmr)

#choosing the response variable 
Y.mdmr = read.csv("y_bt.mdmr.csv",header = T, row.names = 1)
dim(Y.mdmr)

#Transforming data and testing for normality
Y.mdmr_transformed <- log2((Y.mdmr)+1)
apply(Y.mdmr_transformed,2,shapiro.test)
lapply(Y.mdmr_transformed, shapiro.test)
lapply(Y.mdmr_transformed, hist)


#print(Y.mdmr_transformed)
#class(Y.mdmr_transformed)
#str(res)




write.csv(Y.mdmr_transformed, "transformedELISA.csv")
#returns the dimension (e.g. the number of columns and rows) of a matrix, array or data frame
dim(Y.mdmr_transformed)

#All variables are standardized to mean zero and unit variance.

Y.mdmr <- scale(Y.mdmr_transformed,center = TRUE, scale = TRUE)

# Compute distance matrix for the y variable (quantitative/response variable)
D <- dist(Y.mdmr, method = "euclidean")
D

#A Conduct MDMR (regress the quantitative/response variable unto the explanatory variables)
#Option A: Passing the Distance Matrix - note the perm.p = False function
mdmr.res <- mdmr(X = X.mdmr, D = D)


# Check results
#summary(mdmr.res) (You can do this too but this doesn't show results as data frame)

sum <- as.data.frame(summary(mdmr.res))

write.csv(sum, "summary_results.csv")

#Accounting for grouping variable

tempframe <- model.frame(~ Location*Infection_Status*Genotype + Block, data = X.mdmr)
temp.mdmr <- model.matrix(~ Location*Infection_Status*Genotype + Block, data = X.mdmr)
temp.mdmr
write.csv(temp.mdmr, "temp.csv")
res <- mdmr(X = temp.mdmr, D = D)

sum_interaction <- as.data.frame(summary(res))
write.csv(sum_interaction, "summary_interaction_April2024.csv")


