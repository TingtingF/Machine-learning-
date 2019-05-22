# Course material for "Machine Learning in Medical Bioinformatics", 
# Linköping, May 21, 2019.
# Computer Lab (with R): A Cancer Modeling Example.
# See 3_StatPrinciples_RLab.pdf for some background info.

# Starting point for exercise on analysis of miRNA, mRNA and 
# protein data from the paper Aure et al, Integrated analysis 
# reveals microRNA networks coordinately expressed with key 
# proteins in breast cancer, Genome Medicine, 2015.

setwd("C:/Users/tingtinf/Desktop/R")

# Read the data
# sep="\t" tells R that the file is tab-delimited 
# (use " " for space delimited and "," for comma delimited; use "," for a .csv file).

mir = read.table("miRNA-421x282.txt",header = T, sep="\t", dec=".")
rna = read.table("mRNA-100x282.txt", header=T, sep="\t", dec=".")
prt = read.table("prot-100x282.txt", header=T, sep="\t", dec=".")

# Convert to matrix format

mir = as.matrix(mir)
rna = as.matrix(rna)
prt = as.matrix(prt)

# structure of dataset
str(mir)
str(rna)
str(prt)

# Have a look at the data to see if it looks ok

mir[1:4, 1:4]
rna[1:4, 1:4]
prt[1:4, 1:4]

# Look at the overall distribution of expression values

# With the par( ) function, we can include the option mfrow=c(nrows, ncols)
# to create a matrix of nrows x ncols plots that are filled in by row. 
par(mfrow=c(2,2))
# nclass: how many bar should I draw
hist(mir, nclass=40, xlim=c(-5,5), col="lightblue")
hist(rna, nclass=40, xlim=c(-5,5), col="lightblue")
hist(prt, nclass=40, xlim=c(-5,5), col="lightblue")

# Look at mRNA-protein associations (only first nine genes

# The argument pch, an abbreviation for plot character, 
# is the standard argument to set the character that will be plotted in a number of R functions.

par(mfrow=c(3,3))
par(mar=c(3,3,3,3)) #mar – A numeric vector of length 4, which sets the margin sizes in the following order: bottom, left, top, and right. 
for (i in 1:9) {
  plot(rna[i,], prt[i,], pch=".")
}

# Let us add some annotation to the above plot
par(mfrow=c(3,3))
par(mar=c(3,3,3,3))
for (i in 1:9) {
  plot(rna[i,], prt[i,], pch=".")
  title(rownames(rna)[i])
  lines(smooth.spline(rna[i,], prt[i,], df=4), col="red")
}

### TASK: extend the above analysis to cover all genes!
par(mfrow=c(13,8))
par(mar=c(1,1,1,1))
for (i in 1:nrow(rna)) {
  plot(rna[i,], prt[i,], pch=".")
  title(rownames(rna)[i])
  lines(smooth.spline(rna[i,], prt[i,], df=4), col="red")
}
# Compute and plot mRNA-protein correlations

rho = rep(NA, nrow(rna)) # repeat multiple NULL 
for (i in 1:nrow(rna)) {
  rho[i] = cor(rna[i,], prt[i,])
}
par(mfrow=c(1,1))
hist(rho, col="blue")

# Calculate the correlation of each miRNA to each protein

RHO = matrix(NA, nrow(mir), nrow(prt))
for (i in 1:nrow(mir)) {
  for (j in 1:nrow(prt)) {
     RHO[i,j] = cor(mir[i,], prt[j,]) 
  }
}

# Visualize as heatmap

source("clustermap_beta.R") # cluster_map: A master function to perform the full workflow of ClusterMap 

plot.init(tree=c(2,3))
hcluster(RHO, clust="row", distance="euclidean", linkage="complete")
hcluster(RHO, clust="col", distance="euclidean", linkage="complete")
plot.hmap(RHO)
plot.tree(side=2)
plot.tree(side=3)
plot.hmap.key()

#########################################################################

# Model (on the log-scale) the association of miRNA espression on protein 
# expression adjusting for the corresponding mRNA

# Example: Investigate miR-107 and B-RAF (see Aure et al, 2015, Figure 2H):

prt.BRAF = prt[12,]
rna.BRAF = rna[12,]
mir.107 = mir[16,] 

# (a) Linear regression model (on the log-scale) (Aure et al. 2015, equation (3)):
fitA <- lm(prt.BRAF ~ mir.107 + rna.BRAF)
summary(fitA)

#Add smooth non-linear cures to the scatterplots: use existing panel.smooth() function
#Add linear regression lines to the scatterplots:
# abline: to add straight lines to a plot
panel.linear <- function (x, y, col.regres = "blue", ...) 
{ 
  points(x, y) 
  ok <- is.finite(x) & is.finite(y) 
  if (any(ok)) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
} 
pairs(data.frame(mir.107, prt.BRAF, rna.BRAF), 
      lower.panel = panel.smooth,
      upper.panel = panel.linear)


# (b) Lasso-penalised linear model with all miRNAs (Aure et al. 2015, equation (4)):

library(glmnet)

# 10-fold CV to determine the optimal lambda:
# Note: rna.BRAF is penalised together with all the mir variables. Use the penalty.
# factor option to avoid this.

set.seed(1234)

# Glmnet is a package that fits a generalized linear model via penalized maximum likelihood.
# cv.glmnet() uses cross-validation to work out how well each model generalises
# We can automatically find a value for lambda that is optimal by using cv.glmnet()
cvfit <- cv.glmnet(y=prt.BRAF, x=t(rbind(mir, rna.BRAF)),nfolds = 10)

plot(cvfit)

# The lowest point in the curve indicates the optimal lambda: 
# the log value of lambda that best minimised the error in cross-validation. 
# We can extract this values as:
lambda.opt <- cvfit$lambda.min

# Coefficient path plot and coefficients for optimal lambda:
# extract all of the fitted models:
# glmnet.fit a fitted glmnet object for the full data.
fitB <- cvfit$glmnet.fit
summary(fitB)

plot(fitB, xvar="lambda")
abline(v=log(lambda.opt))

# coef is a generic function which extracts model coefficients from objects returned by modeling functions. 

coef(fitB, s=lambda.opt)

# Compare the regression coefficient of mir.107 from the models in (a) and (b):
coef(fitA)["mir.107"]
as.matrix(coef(fitB, s=cvfit$lambda.min))["hsa-miR-107",]


 






