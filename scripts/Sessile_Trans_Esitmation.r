library(phytools)
library(geiger)
library(corHMM)


#Read in data and tree 
phylo.IDs<-read.csv("../data/phylo.IDs.csv",header=TRUE,row.names=1) # classifiers
tree<-read.tree("../data/tree.tre") # pruned tree of 93 taxa
shape<-read.csv("../data/Shape.Data.csv",header=TRUE,row.names=1) # species mean shape data

#Subset Sessile (cementing and nestling) species ("include") from everyone else ("rest")
gp <- rep("rest", length(phylo.IDs$habit))
gp[which(phylo.IDs$habit == "cement" | phylo.IDs$habit == "nestle")] <- "sessile"
gp <- as.factor(gp)

#ML analysis

## Format data
gp_formatted <- cbind(row.names(phylo.IDs),gp)
row.names(gp_formatted) <- row.names(phylo.IDs)

## Plot data with color key for each state
plot.phylo(tree,type = "fan")

colorkey <- c('orange','black')
names(colorkey) <-c("sessile","rest") 
cols<-colorkey[as.character(gp)]

tiplabels(pch = 17,col = cols,cex=0.8)


## Create ER and ARD character trait models
gp_model_er <- getStateMat4Dat(gp_formatted,"ER")
gp_model_ard <- getStateMat4Dat(gp_formatted,"ARD")

#Look at your data 
gp_model_er
gp_model_ard

##plotMKmodel will not work in R 4.2.0
###plotMKmodel(gp_model_er$rate.mat,rate.cat = 1) 
###plotMKmodel(gp_model_ard$rate.mat,rate.cat = 1)

#Fit the models
fit_er <-corHMM(phy = tree, data = gp_formatted, rate.cat = 1, 
                       rate.mat = gp_model_er$rate.mat)
fit_ard <-corHMM(phy = tree, data = gp_formatted, rate.cat = 1, 
                        rate.mat = gp_model_ard$rate.mat)

##The estimated rates
fit_er$solution
fit_ard$solution

#compare models: logL and AIC
fit_er$loglik
fit_ard$loglik

##We can formally perform a likelihood ratio test of nested models by doing the following
LRT <- -2 * (fit_er$loglik - fit_ard$loglik) ## -2 * (model_constrained - model_full)
##p-value
pchisq(LRT,df = 1,lower.tail=FALSE ) ##degrees of freedom is the difference in the number of parameters between the two models. In our case 2-1=1

##Perform ancestral character estimation

##Format the groups
gp.anc <-phylo.IDs$habit
gp.anc[which(phylo.IDs$habit == "cement" | phylo.IDs$habit == "nestle")] <- "sessile"
gp.anc <- as.factor(gp.anc)

anc.ML<-ace(x = gp, phy = tree,  
            type = "discrete", method = 'ML',
            model = gp_model_er$rate.mat)

#Plot the anc
colorkey <- c("#000000", "#E69F00") ##A colorblind friendly color palette
names(colorkey) <- colnames(anc.ML$lik.anc) ##assign each character state a color
tip_cols<-colorkey[as.character(gp.anc)] ##putting the trait values in as.character() is important!

plot(tree) #,type="fan")
nodelabels(pie=anc.ML$lik.anc,piecol=colorkey,cex=0.5)
tiplabels(pch=19,col=tip_cols)
legend(x='bottomleft',legend = names(colorkey),fill=colorkey)





#Bayesian Analysis - run on Terminal
write.nexus(tree, file = "../data/tree.nex")

write.table(gp_formatted, quote = F, row.names = F, col.names = F, 
            file = "../data/sessile.txt")
