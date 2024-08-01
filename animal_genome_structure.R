########## R script used to produce plots for
########## Lewin, Liao and Luo 2024
########## Conservation of animal macrosynteny is the exception not the rule

#### Initially created: 17/07/2024 (Thomas D. Lewin)
#### Last edited: 01/08/2024 (Thomas D. Lewin)

################################# Load packages ################################

library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(ape)
library(phytools)

############################### Load main dataset ##############################

# Load in Dataset S1
ri <- data.frame(read.table("Dataset_S1.txt", sep="\t", header = TRUE))

################################ Fig. 1a,b - PCA ###############################

# Load in dataset with individual SI and CI for each ALG for each species
scores <- read.table("data_for_R_PCA.txt", sep="\t", header = TRUE)

# Make gene column into row name
rownames(scores) <- scores[, 1]
scores$X <- NULL

# transpose the dataframe
scores <- data.frame(scores)
scores_t <- data.frame(t(scores))

# Normalise the dataframe
scores_norm <- scale(scores_t)

# Run PCA
res.pca <- prcomp(scores_norm)

# Plot PCA
fviz_pca_ind(res.pca)
fviz_pca_biplot(res.pca)

# Get contributions of variables to PC1 and PC2
var <- get_pca_var(res.pca)
fviz_contrib(res.pca, choice = "var", axes = 1, top = 20)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 20)

######################## Fig. 1c - High vs Low indicies ########################

## Create boxplot

# Splitting Index 
ggplot(ri, aes(x=Type, y = SI))+ theme_classic() + geom_boxplot()+ylim(0,1)

# Combining Index
ggplot(ri, aes(x=Type, y = CI))+ theme_classic() + geom_boxplot() +ylim(0,1)

## Perform t-tests

# Filter the data for SI
high_group <- ri[ri$Type == "High", "SI"]
low_group <- ri[ri$Type == "Low", "SI"]

# Perform the t-test
t_test_result <- t.test(high_group, low_group)

# get exact p-value
t_test_result$p.value

# Filter the data for CI
high_group <- ri[ri$Type == "High", "CI"]
low_group <- ri[ri$Type == "Low", "CI"]

# Perform the t-test
t_test_result <- t.test(high_group, low_group)

# get exact p-value
t_test_result$p.value

######################### Fig. 1d - RI for each phylum ######################### 

# Reorder species
order <- c("Xenacoelomorpha", "Echinodermata", "Hemichordata", "Chordata",
           "Nematoda", "Nematomorpha", "Tardigrada", "Arthropoda", "Bryozoa",
           "Platyhelminthes", "Annelida", "Mollusca", "Nemertea", "Brachiopoda",
           "Rotifera")

# Set the order of Phylum
ri$Phylum <- factor(ri$Phylum, levels = order)

# Create boxplot 
ggplot(ri, aes(x=Phylum, y = RI, color = Phylum))+ theme_classic() + geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black")

################### Fig. 1e - Ancestral state reconstruction ###################

# This script uses an MCMC approach to ASR called stochastic character mapping
# from Huelsenbeck at al. (2003).

# read in tree
tree <-read.tree("tree.nwk")

# read in file with RI states
X <-read.csv("states_for_ASR.txt",row.names=1, sep = "\t")

# Select column with the state you want to reconstruct - change the number to select
x<-setNames(X[,4],rownames(X))

# Set colours
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
sort(tree$tip.label)
sort(names(x))
state_colours <- c("Low" = "#80B1D3", "High" = "#FB8072")
unique_states <- sort(unique(x))
cols <- state_colours[unique_states]

# Set transition matrix
mod=matrix(c(0,1,0,0),2)

# Make and plot a single simulation
mtree<-make.simmap(tree,x,model=mod)
plotSimmap(mtree,cols,fsize=0.8,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

# A single stochastic character map in isolation is pretty meaningless.
# A whole distribution from a sample of stochastic maps is required. 

# Generate 100 stochastic character maps
mtrees<-make.simmap(tree,x,model=mod,nsim=100)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")

### summarise the maps
pd<-describe.simmap(mtrees,plot=FALSE)

### produce the overall plot
par(mfrow=c(1,1))
plot(pd,fsize=0.6,ftype="i", colors=cols)

########################## Fig. 2a - Oxford dot plots ########################## 

# Produced with SyntenyFinder: see https://github.com/symgenoevolab/SyntenyFinder

######################### Fig. 2b - RI vs seq evolution ######################## 

# Plot the data
plot(ri$BranchLength, ri$RI, xlab="Branch Length", ylab="RI", pch = 19)

## create logistic model
model_logistic <- nls(RI ~ SSlogis(BranchLength, Asym, xmid, scal), data=ri)

# Generate predicted values
predicted_values_logistic <- predict(model_logistic, newdata=data.frame(BranchLength=sort(ri$BranchLength)))

# Plot the fitted logistic curve
lines(sort(ri$BranchLength), predicted_values_logistic, col="black")

########################### Fig. 2c - Idiogram plots ########################### 

# Produced with SyntenyFinder: see https://github.com/symgenoevolab/SyntenyFinder

################################ End of script ################################# 
