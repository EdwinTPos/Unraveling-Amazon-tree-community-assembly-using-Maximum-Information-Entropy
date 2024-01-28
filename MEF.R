#Put the working directory here and clear memory
rm(list=ls())
setwd("yourdirectory")

library(vegan)
library(FD)
library(mice)
library(ggplot2)
library(viridis)
library(doParallel)
library(reshape2)
library(fields)
library(itertools)
library(dismo)

source("MEF_functions.R")	#Function MEE

###########################
##### Field data input ####
###########################

data(BCI)
data(BCI.env)
vegtab = BCI
data_env = BCI.env
traits = read.csv("BCI_trait_data.csv", header = TRUE, row.names = 1)
traits[traits == -99] <- NA
traits = traits[,c(4:10)]

###############################
##### END Field data input ####
###############################



#########################
#### MICE IMPUTATION ####
#########################

traits_imputed = mice(traits,m=5,maxit=10,meth='pmm',seed=500)
summary(traits_imputed)
completedData <- complete(traits_imputed,1)
rownames(completedData) = rownames(traits)
# write.csv(completedData, "BCI_trait_data_imp.csv") #write new .csv file with imputed trait data for easy loading in subsequent runs or analyses.


#############################################################################################
#### Testing the imputations by using the delta approach by Gerko Vink & Stef van Buuren ####
#############################################################################################

#First do a dry run
ini <- mice(traits, maxit = 0)
ini$nmis

#Inspect the missing data
md.pattern(traits)
fx <- fluxplot(traits, cex=.7); fx

#Choose the delta wisely, so look at the summary statistics of the variable to change
summary(na.omit(traits$Sun_LaminaTough)) #what makes up almost a quarter of the average?
delta <- c(0, -100, -200, -300, -400)
imp.all <- vector("list", length(delta))
post <- ini$post

for (i in 1:length(delta)){
  d <- delta[i]
  cmd <- paste("imp[[j]][,i] <- imp[[j]][,i] +", d)
  post["Sun_LaminaTough"] <- cmd
  imp <- mice(traits, post = post, maxit = 15, seed = i, print = FALSE)
  imp.all[[i]] <- imp}

#inspect the imputations for the different delta's
bwplot(imp.all[[1]])
#or do it by a density plot

pdf("mice_delta_adjustment.pdf", paper="a4")
densityplot(imp.all[[1]], lwd = 3)
densityplot(imp.all[[2]], lwd = 3)
densityplot(imp.all[[3]], lwd = 3)
densityplot(imp.all[[4]], lwd = 3)
densityplot(imp.all[[5]], lwd = 3)
dev.off()

## Plot original and estimated data from deviations
plot(complete(imp.all[[1]])$Sun_LaminaTough,complete(imp.all[[1]])$LEAFTHCK_AVD, col="blue", pch=4, xlab="Sun_LaminaTough", ylab = "LEAFTHCK_AVD")
fit.lm = lm(complete(imp.all[[1]])$LEAFTHCK_AVD~complete(imp.all[[1]])$Sun_LaminaTough)
abline(fit.lm, col="blue")
summary(fit.lm)

points(complete(imp.all[[1]])$Sun_LaminaTough,complete(imp.all[[1]])$LEAFTHCK_AVD, col="red", pch=4)
fit.lm2 = lm(complete(imp.all[[5]])$LEAFTHCK_AVD~complete(imp.all[[5]])$Sun_LaminaTough)
abline(fit.lm2, col="red")
summary(fit.lm2)

points(traits$Sun_LaminaTough,traits$LEAFTHCK_AVD, pch=16, ylab = "LEAFTHCK_AVD", xlab="Sun_LaminaTough")
fit.lm3 = lm(traits$LEAFTHCK_AVD~traits$Sun_LaminaTough)
abline(fit.lm3)
summary(fit.lm3)

legend("topleft",pch=c(4, 4, 16), col=c("blue","red","black"), c("Imputed delta = 0","Imputed delta = -400", "Observed"), cex=.8)

#Test the imputations by performing a linear regression between N and P
store.lm = list()
for(i in 1:5){
Sun_LaminaTough_imp = complete(imp.all[[i]])$Sun_LaminaTough
LEAFTHCK_AVD_imp = complete(imp.all[[i]])$LEAFTHCK_AVD
store.lm[[i]] = lm(LEAFTHCK_AVD_imp~Sun_LaminaTough_imp)}

regression_summary = matrix(nrow=7,ncol=7)
for (i in 1:length(store.lm)){
x = cbind(
summary(store.lm[[i]])$coefficients[1,1], #intercept
summary(store.lm[[i]])$coefficients[1,2], #std.error
summary(store.lm[[i]])$coefficients[1,3], #t value
summary(store.lm[[i]])$coefficients[1,4], #Pr
summary(store.lm[[i]])$r.squared,
summary(store.lm[[i]])$adj.r.squared,
summary(store.lm[[i]])$sigma)
regression_summary[i,] = x}

regression_summary[6,]=apply(na.omit(regression_summary), 2, mean)
regression_summary[7,]=apply(na.omit(regression_summary), 2, sd)

rownames(regression_summary) = c(paste("delta",delta), "mean","sd.dev")
colnames(regression_summary) = c("intercept","std.error","t value","Pr","r.squared","adj.r.squared","sigma")
regression_summary


##################################################################################
###### RUNNING THE MAXENT AND MAKE SOME NICE FIGURES (metapop is all plots) ######
##################################################################################

#Run parallel script, clustertype = "FORK" for unix systems and "PSOCK" for windows
run.maxent.par(data=vegtab,traitdata="BCI_trait_data_imp", perm=3, clustertype="FORK")
# save(total.output, file = "output_maxent_alltraitsmice.RData") #save the results from MEF for easy loading in subsequent runs or analyses.

##Load results
load("output_maxent_alltraitsmice.RData") # "output_maxent_alltraitsmice.RData" is the file generated by the function run.maxent.par containing all results


###################################################
#### Merge output from parallelization process ####
###################################################

#Merge all results from the parallel output file
merge.results <<- data.frame()
for(i in 1:length(total.output)){
if (i == 1){merge.results <<- data.frame(total.output[[i]][[1]])}
	else {merge.results<<-rbind(merge.results,data.frame(total.output[[i]][[1]]))}}

#Take out one of the total observed abundances, setting for either [[1]] or [[2]] is the same
merge.ra = total.output[[1]]$observed.ra

#Merge results for predicted ra uniform prior + traits
merge.uniform.ra = data.frame()
for(i in 1:length(total.output)){
if (i == 1){merge.uniform.ra = data.frame(total.output[[i]]$predicted.ra.uniform.plus.traits)}
	else {merge.uniform.ra=rbind(merge.uniform.ra,data.frame(total.output[[i]]$predicted.ra.uniform.plus.traits))}}

#Merge results for predicted ra observed prior + traits
merge.prior.ra = data.frame()
for(i in 1:length(total.output)){
if (i == 1){merge.prior.ra = data.frame(total.output[[i]]$predicted.ra.prior.plus.traits)}
	else {merge.prior.ra=rbind(merge.prior.ra,data.frame(total.output[[i]]$predicted.ra.prior.plus.traits))}}


##############################################
#### Creating Rank Abundance Distribution ####
##############################################

#Fielddata
n.ind.spec	=	colSums(merge.ra)
ME.RAD		=	rev(sort(n.ind.spec))
ME.RAD.rel	=	ME.RAD/sum(ME.RAD)
ranks		  =	seq(1:length(ME.RAD.rel))
plot(ranks, ME.RAD.rel,main = "Relative abundance distribution", xlab="rank", ylab="log(proportional abundance)", col =  "green", cex = 0, type="o", xlim=c(0,length(ranks)+10), log="y")

#Only traits fit
n.ind.spec	=	colSums(merge.uniform.ra)
ME.RAD		=	rev(sort(n.ind.spec))
ME.RAD.rel	=	ME.RAD/sum(ME.RAD)
ranks		  =	seq(1:length(ME.RAD.rel))
lines(ranks,ME.RAD.rel,col="blue")

#Traits+Metacommunity fit
n.ind.spec	=	colSums(merge.prior.ra)
ME.RAD		=	rev(sort(n.ind.spec))
ME.RAD.rel	=	ME.RAD/sum(ME.RAD)
ranks		  =	seq(1:length(ME.RAD.rel))
lines(ranks,ME.RAD.rel,col="red")

####################################################################################
#### Plot observed versus predicted relative abundances per plot for each model ####
####################################################################################

tiff(file="ObservedvsPredictedRA.tiff",width = 6, height = 6, units = "in", res=400)
par(mfrow=c(2,2))
##Linear regression test uniform prior + traits
test = as.vector(merge.ra); test2 = as.vector(as.matrix(merge.uniform.ra))
plot(test,test2,xlim=c(0.00005,1),ylim=c(0.00000001,1),log="xy", type="n", xlab="log(observed relative abundance)", ylab="log(predicted relative abundance)", main="uniform+traits per plot")
points(test,test2, col=ifelse(test >= .1,"red","black"))
fit = cor(test,test2,method="pearson"); abline(0,1,untf=T)
legend("bottomright", legend = paste("Pearson R2 = ",round(fit, digits=2)), cex=.7, bty="n")

#Uniform prior plus traits, summed over all plots
plot(colSums(merge.ra),colSums(merge.uniform.ra), log="xy", col=ifelse(colSums(merge.ra) >= 1,"red","black"), xlab="log(observed regional relative abundance)", ylab="log(predicted regional relative abundance)",main="uniform+traits summed", ylim=c(0.0001,100),xlim=c(0.0001,100))
fit = cor(colSums(merge.ra),colSums(merge.uniform.ra),method="pearson"); abline(0,1,untf=T)
legend("bottomright", legend = paste("Pearson R2 = ",round(fit, digits=2)), cex=.7, bty="n")

##Linear regression test observed vs observed prior plus traits
test = as.vector(merge.ra); test2 = as.vector(as.matrix(merge.prior.ra))
plot(test,test2,xlim=c(0.00005,1),ylim=c(0.000000001,1),log="xy", type="n",xlab="log(observed relative abundance)", ylab="log(predicted relative abundance)",main="prior+traits per plot")
points(test,test2, col=ifelse(merge.ra >= .1,"red","black"))
fit = cor(test,test2,method="pearson"); abline(0,1,untf=T)
legend("bottomright", legend = paste("Pearson R2 = ",round(fit, digits=2)), cex=.7, bty="n")

#Observed prior plus traits, summed over all plots
plot(colSums(merge.ra),colSums(merge.prior.ra), log="xy", col=ifelse(colSums(merge.ra) >= 1,"red","black"),xlab="log(observed regional relative abundance)", ylab="log(predicted regional relative abundance)", main="prior+traits summed",ylim=c(0.0001,100),xlim=c(0.0001,100))
fit = lm(colSums(merge.prior.ra)~colSums(merge.ra))
fit = cor(colSums(merge.ra),colSums(merge.prior.ra),method="pearson"); abline(0,1,untf=F)
legend("bottomright", legend = paste("Pearson R2 = ",round(fit, digits=3)), cex=.7, bty="n")

dev.off()
par(mfrow=c(1,1)) #resetting window device settings 

############################################################################################
#### Boxplots concerning relative proportion of total biologically relevant information ####
############################################################################################
results = data.frame(merge.results[,c("pure.trait.effect","pure.metacommunity.effect", "joint.trait.metacommunity.effect","unexplained.effect")])
results$Habitat = 0

#Note that habitats here are used for categories to compare for example, with habitat derived from the environmental data on BCI available from vegan.
for(i in 1:length(unique(data_env$Habitat))){
type = as.character(unique(data_env$Habitat)[i])
selecthabitat = data_env[data_env$Habitat==type,]
results$Habitat[rownames(results) %in% rownames(selecthabitat)] = type}

#Note that the below only loads results with positive values. The MaxEnt functions generates negative KLR2 values when model bias is higher than either pure effects which is not corrected for in the original functions. 
test = results[rowSums(results[,1:4] >= 0)==4,]
test.m <- melt(test)

myplot = ggplot(test.m, aes(x = variable, y = value, fill = Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("blue", "yellow","green","brown","purple", "orange"))

myplot + theme_bw() + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank())

#############################################################
#### Boxplots concerning lambda values per type/habitat  ####
#############################################################

merge.lambdas <<- data.frame()
for(i in 1:length(total.output)){
if (i == 1){merge.lambdas <<- data.frame(total.output[[i]][[2]])}
	else {merge.lambdas<<-rbind(merge.lambdas,data.frame(total.output[[i]][[2]]))}}

colnames(merge.lambdas)
lambda.select=merge.lambdas[2:8]
lambda.select = scale(lambda.select)
results = data.frame(lambda.select)
results$Habitat = 0

for(i in 1:length(unique(data_env$Habitat))){
type = as.character(unique(data_env$Habitat)[i])
selecthabitat = data_env[data_env$Habitat==type,]
results$Habitat[rownames(results) %in% rownames(selecthabitat)] = type}

results= results[!results$Habitat == "Young",]
results= results[!results$Habitat == "Swamp",]
data = melt(results)

anova_table=data.frame(matrix(nrow=length(unique(data$variable)), ncol=2))
x = c("constraint","label")
colnames(anova_table)=x

for(i in 1:length(unique(data$variable))){
anova_table[i,"constraint"] = as.character(unique(data$variable)[i])
anova = aov(value~Habitat, data=subset(data, variable == unique(data$variable)[i]))
if(summary(anova)[[1]]$P[1] <= .05) stars="*"
if(summary(anova)[[1]]$P[1] <= .001) stars="**"
if(summary(anova)[[1]]$P[1] <= .0001) stars="***"
label = paste("F=",round(summary(anova)[[1]]$F[1],3),stars, sep="")
anova_table[i,"label"] =label }

##Barplots
agg = aggregate(results, by =list(results$Habitat), FUN=mean)
test_agg = agg[,1:8]
data = melt(test_agg)

myplot = ggplot(data, aes(fill=Group.1, y=value, x=variable)) +
    geom_bar(position="dodge", stat="identity")+
	
  	scale_fill_manual(values = c("blue","yellow2","green","brown","purple","orange"), name = "Habitat")+
	labs(x=" ",y="Lambda value")+
	annotate("text", x=1, y=2.5, label= anova_table$label[1], cex=3)+
	annotate("text", x=2, y=2.5, label= anova_table$label[2], cex=3)+
	annotate("text", x=3, y=2.5, label= anova_table$label[3], cex=3)+
	annotate("text", x=4, y=2.5, label= anova_table$label[4], cex=3)+
	annotate("text", x=5, y=2.5, label= anova_table$label[5], cex=3)+
	annotate("text", x=6, y=2.5, label= anova_table$label[6], cex=3)+
	annotate("text", x=7, y=2.5, label= anova_table$label[7], cex=3)


myplot + theme_bw() + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.x=element_blank(),axis.title.x=element_blank())+geom_vline(xintercept=seq(1.5, by=1, 14), linetype="dashed", color="blue", size=.2)

#####################################
#### Create table for all values ####
#####################################
data.output = data.frame(merge.results)
colnames(data.output)

forest.types = data_env$Habitat
unique(forest.types)

fulltable = matrix(nrow=14, ncol = length(unique(forest.types))+1)
rownames(fulltable) = colnames(data.output)
colnames(fulltable) = c(as.character(unique(forest.types)),"Combined")

for(i in 1:length(unique(forest.types))){
type = as.character(unique(forest.types)[i])
selectforesttype = data_env[forest.types==type,]
data.output2 = data.output[rownames(data.output) %in% rownames(selectforesttype),]
data.output3 = apply(data.output2,2,mean)
fulltable[,i] = data.output3
}

fulltable[,6] = apply(data.output,2,mean)
rownames(fulltable) = c("Rkl(u)","Rkl(u,t)","Rkl(m)","Rkl(m,t)","ΛR2KL(t|φ)","ΛR2KL(m|φ)","ΛR2KL(t|m)","ΛR2KL(m|t)","ΛR2KL(m+t)","1-ΛR2KL(m+t)","Pure trait effect","Pure metacommunity effect","Joint effect","Unexplained effects")

#Generate table with all results
fulltable

#########################################################################
##### Running the MEF procedure for incrementing metacommunity size #####
#########################################################################

vegtab.field = vegtab
utms <- SpatialPoints(data_env[, c("UTM.EW", "UTM.NS")],
                      proj4string=CRS("+proj=utm +zone=16"))
longlats <- spTransform(utms, CRS("+proj=longlat"))
longlat = data.frame(longlats); colnames(longlat) = c("lon","lat")

##remove empty columns
vegtab.field = vegtab.field[,colSums(vegtab.field) != 0]
n.spec = length(colnames(vegtab.field)); n.spec
length(rownames(vegtab.field)) == length(data_env$Habitat) #check if all plots are in

##Calculate geographical distances and base range of categories on results
p <- data.frame(lat=longlat$lat, lon=longlat$lon)
rough_distances <<- rdist.earth(p, miles=F)
range(rough_distances)
categories = seq(0.01,1, by=0.1) #set categories based on range of geographical distances
categories = seq(0.01,1, by=0.25) #set categories based on range of geographical distances

## Run the actual looping function to increase metacommunity size subsequently in parallel, again use clustertype = "FORK" for unix systems and "PSOCK" for windows 
time = proc.time()
meta.estimate.par(vegtab.field,longlat$lat,longlat$lon,traitdata="BCI_trait_data_imp",
						categories,perm=2, clustertype="FORK", cutoff=.10)
proc.time()-time


######################################################################################################
######################################################################################################
#NOTE: function meta.estimate.par() makes use of a connection:
# conn <- file(sprintf(""filepath"/output_%d.txt" , Sys.getpid()) , open = "a" )

# after the calculations have been done, the files need to be merged and purged for duplicates using a terminal and working in the specific folder:
#for f in *.txt ; do sort $f | uniq >> short_data.txt ;done 

######################################################################################################
######################################################################################################

###########################
### LOADING THE RESULTS ###
###########################

#Load already trimmed data
effects_distance = read.csv("short_data.csv",header=F, encoding="latin1")
colnames(effects_distance)= c("Plotcode", "maxdistance","pure.trait.effect","pure.metacommunity.effect", "joined.effect","unexplained.effect")
effects_distance = na.omit(effects_distance)
#Number of plots
length(unique(effects_distance$Plotcode))

###########################################################
###Fit non-linear regression and use exponent for maps ####
###########################################################

#Use all points and determine exponent from simple exponential regression
metatrait_dist_ratio = effects_distance[,c(1,2,4)]
ratio_metacom = matrix(, ncol=2, nrow = length(unique(metatrait_dist_ratio[,1])))

	for (i in 1:length(unique(metatrait_dist_ratio$Plotcode))){
		name = unique(metatrait_dist_ratio$Plotcode)[i]
		y = metatrait_dist_ratio[metatrait_dist_ratio$Plotcode == name,3]
		x = metatrait_dist_ratio[metatrait_dist_ratio$Plotcode == name,2]
		df = data.frame(x,y); df = df[order(x),]
		df = as.matrix(df)
		
		fit_nls = nls(df[,2] ~ I(df[,1]^power), start = list(power = 1),	
		control=nls.control(maxiter=10000))

		#z = predict(fit_nls,categories)
		#ratio = z[49]/z[1]
		ratio = summary(fit_nls)$parameters[1]
		ratio_metacom[i,1] = name
		ratio_metacom[i,2] = ratio
			}

ratio_metacom = data.frame(ratio_metacom)
colnames(ratio_metacom) = c("PlotCode","Exp")
ratio_metacom$PlotCode = as.character(unique(metatrait_dist_ratio$Plotcode))
ratio_metacom = cbind(ratio_metacom,data_env[rownames(data_env) %in%ratio_metacom$PlotCode, c("UTM.EW","UTM.NS")])

ggplot(ratio_metacom) + 
geom_tile(aes(x=UTM.EW, y = UTM.NS, fill=Exp*1000)) + 
  coord_fixed(ratio = 1) +
  scale_fill_viridis(direction = -1)


##############################################################
### Plot distace category versus pure metacommunity effect ###
##############################################################

#Note that in this version there is no categorization done for forestplots given the small spatial scale for BCI but all plots are plotted overall.

##Plot loess prediction based on all points
test		= cbind(effects_distance$maxdistance, effects_distance$pure.metacommunity.effect)
test2	= test[order(effects_distance$maxdistance),]
test3	= aggregate(test2,list(test2[,1]),FUN="mean")

a <- which(agg_selection$maximum==min(agg_selection$maximum),arr.ind=T)
b <- effects_distance$Plotcode == agg_selection[a,]$names.i.
min.effect = effects_distance$pure.metacommunity.effect[b]
min.effect.distance = effects_distance$maxdistance[b]

a <- which(agg_selection$difference==max(agg_selection$difference),arr.ind=T)
b <- effects_distance$Plotcode == agg_selection[a,]$names.i.
max.effect = effects_distance$pure.metacommunity.effect[b]
max.effect.distance = effects_distance$maxdistance[b]

#LOESS prediction MIN
myloess_min = loess(min.effect~min.effect.distance, na.action = na.exclude, span=.10) 
prediction_min = predict(myloess_min, min.effect.distance, se=T)

#LOESS prediction MAX
myloess_max = loess(max.effect~max.effect.distance, na.action = na.exclude, span=.10) 
prediction_max = predict(myloess_max, max.effect.distance, se=T)

#LOESS prediction
myloess <- loess(effects_distance$pure.metacommunity.effect~effects_distance$maxdistance, na.action = na.exclude, span=.4) 
prediction = predict(myloess, test3$V1, se=T)

#Create plot with polygon based on min and max of prediction SE
plot(test3$V1,test3$V2, xlab="Distance from focal plot (km)", ylab="Pure metacommunity effect")
lines(test3$V1, (prediction$fit-prediction_min$se), lty="dashed", col="lightskyblue1", lwd=1) 
lines(test3$V1, (prediction$fit+prediction_min$se), lty="dashed", col="lightskyblue1", lwd=1) 
polygon(c(test3$V1, rev(test3$V1)), c((prediction$fit-prediction_min$se), rev((prediction$fit+prediction_min$se))), col="lightskyblue1", border=NA)
points(test3$V1,test3$V2, col="black", pch=20)
lines(test3$V1,prediction$fit, col = "red", lwd = 1.5)


######################################################
### Plot distace category versus pure trait effect ###
######################################################

#Note that this does include the negative values for the trait values, as indicated above an artefact from the model bias giving more information than either pure effects. 

test = cbind(effects_distance$maxdistance, effects_distance$pure.trait.effect)
test2 = test[order(effects_distance$maxdistance),]
test3=aggregate(test2,list(test2[,1]),FUN="mean")

a <- which(agg_selection$maximum==min(agg_selection$maximum),arr.ind=T)
b <- effects_distance$Plotcode == agg_selection[a,]$names.i.
min.effect = effects_distance$pure.trait.effect[b]
min.effect.distance = effects_distance$maxdistance[b]

a <- which(agg_selection$difference==max(agg_selection$difference),arr.ind=T)
b <- effects_distance$Plotcode == agg_selection[a,]$names.i.
max.effect = effects_distance$pure.trait.effect[b]
max.effect.distance = effects_distance$maxdistance[b]

#LOESS prediction MIN
myloess_min = loess(min.effect~min.effect.distance, na.action = na.exclude, span=.10) 
prediction_min = predict(myloess_min, min.effect.distance, se=T)

#LOESS prediction MAX
myloess_max = loess(max.effect~max.effect.distance, na.action = na.exclude, span=.10) 
prediction_max = predict(myloess_max, max.effect.distance, se=T)

#LOESS prediction
myloess <- loess(effects_distance$pure.trait.effect~effects_distance$maxdistance, na.action = na.exclude, span=.4) 
prediction = predict(myloess, test3$V1, se=T)

#Create plot with polygon based on min and max of prediction SE
plot(test3$V1,test3$V2, xlab="Distance from focal plot (km)", ylab="Pure trait effect")
axis(side = 1, at = categories)
lines(test3$V1, (prediction$fit-prediction_min$se), lty="dashed", col="lightskyblue1", lwd=1) 
lines(test3$V1, (prediction$fit+prediction_min$se), lty="dashed", col="lightskyblue1", lwd=1) 
polygon(c(test3$V1, rev(test3$V1)), c((prediction$fit-prediction_min$se), rev((prediction$fit+prediction_min$se))), col="lightskyblue1", border=NA)
points(test3$V1,test3$V2, col="black", pch=20)
lines(test3$V1,prediction$fit, col = "red", lwd = 1.5)
