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
D <- dist(Y.mdmr, method = "manhattan")
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
write.csv(sum_interaction, "summary_interaction.csv")

#mixed.res <- mixed.mdmr(~ Genotype * Location * Infection_Status * Block + (Genotype + Location + Infection_Status ||Genotype), data = X.mdmr, D = D)

#sum2 <- as.data.frame(summary(mixed.res))
#write.csv(sum2, "summary_results2.csv")

# Study univariate effect sizes
#The first option to compute ?? statistics it to pass the matrix of predictors, 
#matrix of outcome data, and distance metric to the function delta(). The additional argument niter 
#determines how many times the pseudo-jackknife procedure is conducted and averaged, such that larger 
#numbers yield more stable estimates but require more computation time)
delta.res <- delta(temp.mdmr, Y = Y.mdmr, dtype = "euclidean", 
                   niter = 5, plot.res = T)


#Option B: Passing Gower’s Trannsformed Distance Matrix
#G <- gower(D)
#lambda <- eigen(G, only.values = T)$values
#mdmr.res2 <- mdmr(X = X.mdmr, G = G, lambda = lambda)

#sum <- as.data.frame(summary(mdmr.res2))

#D <- dist(Y.mdmr, method = "manhattan")
#G <- gower(D)
#q <- ncol(Y.mdmr)
#G.list <- lapply(1:q, FUN = function(k) {
 # set.seed(k)
  #Y.shuf <- Y.mdmr
  #Y.shuf[,k] <- sample(Y.shuf[,k])
  #gower(dist(Y.shuf, method = "manhattan"))
#})
#names(G.list) <- colnames(Y.mdmr)
#par(mar = c(5, 5, 4, 2) + 0.1)
#delta(X = X.mdmr, G = G, G.list = G.list, plot.res = T)
