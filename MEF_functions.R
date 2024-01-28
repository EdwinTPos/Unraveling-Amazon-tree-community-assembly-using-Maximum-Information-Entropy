##########################################################################
###### FUNCTION FOR RUNNING MAXENT PROCEDURE FOR EACH PLOT PARALLEL ######
##########################################################################

run.maxent.par = function(data, traitdata, perm, clustertype){

Nperm = perm

#Reading in the total data
abund = data
abund = abund[,colSums(abund) != 0]
abund = abund[,sort(colnames(abund))]

#Setting up the traits
traits = read.csv(traitdata,header=T)
rownames(traits) = traits[,1]
traits = traits[,-1]
traits = traits[rownames(traits) %in% colnames(abund), ]
traits = traits[order(rownames(traits)),]

#Setting up the total abundance data, use only species left from traits after removal because of NA
abund <- abund[ ,colnames(abund) %in% rownames(traits)] # abundance matrix
abund <- abund+1e-4
abund <- t(apply(abund, 1, function(x) x / sum(x) )) # relative

#Calculate total community weighted means
CWM <- functcomp(traits, abund) # community-aggregated traits
traits <- t(traits) # transpose matrix
ra = abund

#Set the Prior
q = colSums(abund)/sum(abund)

ncores <- detectCores()
cl <- makeCluster(ncores, type=clustertype)
registerDoParallel(cl)

results = foreach(m=isplitRows(data, chunks=ncores),.packages='FD', .export=c("maxent2", 					"maxent.test2")) %dopar% { 

Nsims = length(rownames(m))
plot.set=F

#Set the different parameters to NA to be filled during the loop
model.mean.null.given.uniform <- model.uniform.prior.plus.traits <- model.mean.null.given.prior<-model.prior.plus.traits<-cor.local.meta<-cor.trait.x<-rep(NA,Nsims)
delta.R.trait<-delta.R.neutral<-rep(NA,Nsims)
information.unique.to.local.trait.constraints<-information.unique.to.neutral.prior<- joint.information<-biologically.unexplained.information<-rep(NA,Nsims)

#Set the matrices to be filled for the trait and local community correlations
cor.meta.x = matrix(0,1,ncol=length(rownames(m)))
rownames(cor.meta.x) = "metacommunity.abundance"
colnames(cor.meta.x) = rownames(m)

cor.trait.x = matrix(0,nrow=length(rownames(traits)),ncol=length(rownames(m)))
rownames(cor.trait.x) = rownames(traits)
colnames(cor.trait.x) = rownames(m)

#Set the matrices to be filled for the different parameters 
model.mean.null.given.uniform = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(model.mean.null.given.uniform) = rownames(m)
colnames(model.mean.null.given.uniform) = "model.mean.null.given.uniform"

model.mean.null.given.prior = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(model.mean.null.given.prior) = rownames(m)
colnames(model.mean.null.given.prior) = "model.mean.null.given.uniform"

model.uniform.prior.plus.traits = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(model.uniform.prior.plus.traits) = rownames(m)
colnames(model.uniform.prior.plus.traits) = "model.uniform.prior.plus.traits"

model.prior.plus.traits = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(model.prior.plus.traits) = rownames(m)
colnames(model.prior.plus.traits) = "model.prior.plus.traits"

lambda.values = matrix(0,nrow = length(rownames(m)),ncol=2+length(rownames(traits)))
rownames(lambda.values) = rownames(m)
colnames(lambda.values) = c("intercept", rownames(traits), "prior")

delta.R.trait = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(delta.R.trait) = rownames(m)
colnames(delta.R.trait) = "delta.R.trait"

delta.R.neutral = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(delta.R.neutral) = rownames(m)
colnames(delta.R.neutral) = "delta.R.neutral"

information.unique.to.local.trait.constraints = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(information.unique.to.local.trait.constraints) = rownames(m)
colnames(information.unique.to.local.trait.constraints) = "information.unique.to.local.trait.constraints"

information.unique.to.neutral.prior = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(information.unique.to.neutral.prior) = rownames(m)
colnames(information.unique.to.neutral.prior) = "information.unique.to.neutral.prior"

joint.information = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(joint.information) = rownames(m)
colnames(joint.information) = "joint.information"

biologically.unexplained.information = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(biologically.unexplained.information) = rownames(m)
colnames(biologically.unexplained.information) = "biologically.unexplained.information"

pure.trait.effect = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(pure.trait.effect) = rownames(m)
colnames(pure.trait.effect) = "pure.trait.effect"

pure.metacommunity.effect = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(pure.metacommunity.effect) = rownames(m)
colnames(pure.metacommunity.effect) = "pure.metacommunity.effect"

joint.trait.metacommunity.effect = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(joint.trait.metacommunity.effect) = rownames(m)
colnames(joint.trait.metacommunity.effect) = "joint.trait.metacommunity.effect"

unexplained.effect = matrix(0,nrow = length(rownames(m)),ncol=1)
rownames(unexplained.effect) = rownames(m)
colnames(unexplained.effect) = "unexplained.effect"

predicted.ra.uniform.plus.traits = matrix(0,nrow = length(rownames(m)),ncol=length(colnames(abund)))
rownames(predicted.ra.uniform.plus.traits) = rownames(m)
colnames(predicted.ra.uniform.plus.traits) = colnames(abund)

predicted.ra.prior.plus.traits = matrix(0,nrow = length(rownames(m)),ncol=length(colnames(abund)))
rownames(predicted.ra.prior.plus.traits) = rownames(m)
colnames(predicted.ra.prior.plus.traits) = colnames(abund)

#Now run the loop for each plot
for(i in 1:Nsims){

#Calculate correlation between local abundances and metacommunity or trait values
cor.local.meta[i]<-cor(as.vector(ra[rownames(m)[i],]),as.vector(q))
cor.meta.x[i] = cor.local.meta[i]

for(j in 1:length(rownames(traits))){
trait.cor<-cor(as.vector(traits[j,]),as.vector(ra[rownames(m)[i],]))
cor.trait.x[j,i] = trait.cor}

#Run first model fit using CWM and only uniform prior
fit1<-maxent2(constr=CWM[rownames(m)[i],],states=traits)
temp1<-maxent.test2(model=fit1,obs=ra[rownames(m)[i],],nperm=Nperm,quick=FALSE, plot=plot.set)
if(plot.set == T)text(0.001,0.001,"uniform prior")
predicted.ra.uniform.plus.traits[i,] = fit1$prob

#An ESTIMATE (using N permutations) of model bias (model.mean.null.given.uniform)
model.mean.null.given.uniform[i,]<-temp1$mean.KLR2.null
model.bias<-model.mean.null.given.uniform[i,]

#Set model fit using the CWM and only the uniform prior
model.uniform.prior.plus.traits[i,]<-temp1$KLR2.prior.plus.traits

#If model using the maximally uninformative prior plus traits is < the model bias, then correct
if(model.uniform.prior.plus.traits[i,]<model.bias)model.uniform.prior.plus.traits[i,]<-model.bias

#Calculate model using the CWM plus neutral prior
fit2<-maxent2(constr=CWM[rownames(m)[i],],states=traits,prior=q,lambda=TRUE)
fitted.lambdas<-fit2$lambda
temp2<-maxent.test2(model=fit2,obs=ra[rownames(m)[i],],nperm=Nperm,quick=FALSE, plot=plot.set)
if(plot.set == T)text(0.001,0.001,"neutral prior")
model.mean.null.given.prior[i,]<-temp2$mean.KLR2.null
predicted.ra.prior.plus.traits[i,] = fit2$prob

#If this null model, given the neutral prior, is < the mean given the maximally uniformative prior, then correct
if(model.mean.null.given.prior[i,]<model.bias)model.mean.null.given.prior[i,]<-model.bias

#Now fit using the CWM and the neutral prior:
model.prior.plus.traits[i,]<-temp2$KLR2.prior.plus.traits

# if this model, given the neutral prior and the traits, is less than the fit of the model using the maximally uniformative prior and the traits, then correct:
if(model.prior.plus.traits[i,]<model.uniform.prior.plus.traits[i,])model.prior.plus.traits[i,]<-model.uniform.prior.plus.traits[i,]

#Now generate the necessary numbers for the output
just.neutral1<-model.prior.plus.traits[i,]-model.uniform.prior.plus.traits[i,]
just.neutral2<-model.mean.null.given.prior[i,]-model.mean.null.given.uniform[i,]
just.traits1<-model.uniform.prior.plus.traits[i,]-model.mean.null.given.uniform[i,]
just.traits2<-model.prior.plus.traits[i,]-model.mean.null.given.prior[i,]

information.unique.to.local.trait.constraints[i,]<-just.traits2
information.unique.to.neutral.prior[i,]<-just.neutral1

delta.R.trait[i,]<-just.traits1
delta.R.neutral[i,]<-just.neutral2
lambda.values[i,]<-fitted.lambdas

T1<-just.traits1-just.traits2
T2<-just.neutral2-just.neutral1

if(abs(T1-T2)>1e-6)return(
 list(message="problemtest",T1=T1,T2=T2,traits1=just.traits1,traits2=just.traits2,neutral1=just.neutral1,neutral2=just.neutral2, diff.traits=abs(just.traits2-just.traits1),diff.neutral=abs(just.neutral1-just.neutral2)))

joint.information[i,]<-T1

biologically.unexplained.information[i,]<-(1-model.prior.plus.traits[i,])

pure.trait.effect[i,]= information.unique.to.local.trait.constraints[i,]/(1-model.mean.null.given.uniform[i,]) 

pure.metacommunity.effect[i,]= information.unique.to.neutral.prior[i,]/(1-model.mean.null.given.uniform[i,]) 

joint.trait.metacommunity.effect[i,]= joint.information[i,]/(1-model.mean.null.given.uniform[i,]) 

unexplained.effect[i,]= biologically.unexplained.information[i,]/(1-model.mean.null.given.uniform[i,]) 


#Finally generate the output itself
output = list(

model.output.per.plot = cbind(model.mean.null.given.uniform, model.uniform.prior.plus.traits, model.mean.null.given.prior, model.prior.plus.traits, delta.R.trait, delta.R.neutral,information.unique.to.local.trait.constraints, information.unique.to.neutral.prior, joint.information,biologically.unexplained.information, pure.trait.effect,pure.metacommunity.effect,joint.trait.metacommunity.effect,unexplained.effect),

lambdas = lambda.values,

model.mean.null.given.uniform=round(c(mean(model.mean.null.given.uniform),
sd(model.mean.null.given.uniform)/sqrt(Nsims)),4),

model.uniform.prior.plus.traits=round(c(mean(model.uniform.prior.plus.traits),
sd(model.uniform.prior.plus.traits)/sqrt(Nsims)),4),

model.mean.null.given.prior=round(c(mean(model.mean.null.given.prior),
sd(model.mean.null.given.prior)/sqrt(Nsims)),4), 

model.prior.plus.traits=round(c(mean(model.prior.plus.traits),
sd(model.prior.plus.traits)/sqrt(Nsims)),4),

mean.pure.trait.effect=round(c(mean(pure.trait.effect),
sd(pure.trait.effect)/sqrt(Nsims)),4),

mean.pure.metacommunity.effect=round(c(mean(pure.metacommunity.effect),
sd(pure.metacommunity.effect)/sqrt(Nsims)),4),

mean.joint.trait.metacommunity.effect=round(c(mean(joint.trait.metacommunity.effect),
sd(joint.trait.metacommunity.effect)/sqrt(Nsims)),4),

mean.unexplained.effect=round(c(mean(unexplained.effect),
sd(unexplained.effect)/sqrt(Nsims)),4),

delta.R.trait=round(c(mean(delta.R.trait),
sd(delta.R.trait)),4),
delta.R.neutral=round(c(mean(delta.R.neutral),
sd(delta.R.neutral)),4),
information.unique.to.local.trait.constraints=round(c(mean(information.unique.to.local.trait.constraints),
sd(information.unique.to.local.trait.constraints)/sqrt(Nsims)),4),

information.unique.to.neutral.prior=round(c(mean(information.unique.to.neutral.prior),
sd(information.unique.to.neutral.prior)/sqrt(Nsims)),4),

joint.information=round(c(mean(joint.information),
sd(joint.information)/sqrt(Nsims)),4),

biologically.unexplained.information=round(c(mean(biologically.unexplained.information),
sd(biologically.unexplained.information)/sqrt(Nsims)),4),

mean.cor.local.metacommunity.ra=round(c(mean(cor.meta.x),
sd(cor.meta.x)/sqrt(Nsims)),4),

cor.local.metacommunity.ra = cor.meta.x, 
mean.cor.trait.x=round(c(mean(cor.trait.x),
sd(cor.trait.x)/sqrt(Nsims)),4),

cor.local.trait.ra = cor.trait.x,

predicted.ra.uniform.plus.traits = predicted.ra.uniform.plus.traits,
predicted.ra.prior.plus.traits = predicted.ra.prior.plus.traits,

observed.ra=ra,observed.metacommunity.ra=q)}
return(output)
}

stopCluster(cl)

total.output <<- results

}






#####################################################################
###### FUNCTION FOR RUNNING DISTANCES FROM FOCAL PLOT PARALLEL ######
#####################################################################


meta.estimate.par = function(data,Lat, Lon,traitdata, categories, perm, clustertype, cutoff){

#Set the dataframes to be filled
p <- data.frame(lat=Lat, lon=Lon)
d <- setNames(do.call(rbind.data.frame,combn(1:nrow(p),2,simplify=F)),c('p1','p2'))

#Calculate geographic distances 
rough_distances <<- rdist.earth(p, miles=F)

#Reading in the total data
total.abund = data
total.abund = total.abund[,colSums(total.abund) != 0]
total.abund = total.abund[,sort(colnames(total.abund))]

#Setting up the traits
traits = read.csv(traitdata,header=T)
rownames(traits) = traits[,1]
traits = traits[,-1]
traits = traits[rownames(traits) %in% colnames(total.abund), ]
traits = traits[order(rownames(traits)),]

#Setting up the total abundance data, use only species left from traits after removal because of NA
total.abund <- total.abund[ ,colnames(total.abund) %in% rownames(traits)] # abundance matrix
total.abund <- total.abund+1e-4
total.abund <- t(apply(total.abund, 1, function(x) x / sum(x) )) # relative

#Calculate community weighted means
CWM.total <- functcomp(traits,total.abund) # community-aggregated traits
traits <- t(traits) # transpose matrix

ncores <- detectCores()
cl <- makeCluster(ncores, type=clustertype)
registerDoParallel(cl)

results = foreach(m=isplitRows(data, chunks=ncores),.packages='FD',.export=c("maxent2", "maxent.test2","rough_distances")) %dopar% { 

conn <- file(sprintf(paste(as.character(getwd()),"/output_%d.txt", sep="") , Sys.getpid()) , open = "a" )

for (k in 1:length(rownames(m))){

stop = FALSE
plot.set=F
Nperm = 2

plotname = rownames(m)[k]
s = which(rownames(data) == plotname)  

#Create subsetted dataframe for plots to compare with
d$distances = 0
d = d[d$p1==s,]
all_plots = seq(1:length(rownames(data)))
d = data.frame(s,all_plots)
colnames(d) = c("p1","p2")

#fill subsetted dataframes, sort by distances
d$distances = rough_distances[,s]
d = d[-s,]
newdata = d[order(d$distances),]
newdata$distances = as.numeric(newdata$distances)

for (i in 1:length(categories)){
range = newdata[newdata$distances<categories[i],]
if(length(range$distances) < 2) next

names1 = unique(range$p1)
names2 = unique(range$p2)
rowstotake = c(names1, names2)

#Reading in the data
abund = total.abund[rowstotake,]
abund = abund[,colSums(abund) != 0]
abund = abund[,sort(colnames(abund))]
ra = abund[1,]

#Take CWM and set prior
CWM <- CWM.total[rownames(abund)[1],]

q = colSums(abund)/sum(abund)

#Set the different parameters to NA to be filled during the loop
model.mean.null.given.uniform <- model.uniform.prior.plus.traits <- model.mean.null.given.prior<-model.prior.plus.traits<-cor.local.meta<-cor.trait.x<-NA
delta.R.trait<-delta.R.neutral<-NA
information.unique.to.local.trait.constraints<-information.unique.to.neutral.prior<- joint.information<-biologically.unexplained.information<-NA

#Set the matrices to be filled for the trait and local community correlations
cor.meta.x = matrix(0,1,ncol=length(rownames(abund)))
rownames(cor.meta.x) = "metacommunity.abundance"
colnames(cor.meta.x) = rownames(abund)

cor.trait.x = matrix(0,nrow=length(rownames(traits)),ncol=length(rownames(abund)))
rownames(cor.trait.x) = rownames(traits)
colnames(cor.trait.x) = rownames(abund)

#Set the matrices to be filled for the different parameters 
model.mean.null.given.uniform = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(model.mean.null.given.uniform) = rownames(abund)
colnames(model.mean.null.given.uniform) = "model.mean.null.given.uniform"

model.mean.null.given.prior = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(model.mean.null.given.prior) = rownames(abund)
colnames(model.mean.null.given.prior) = "model.mean.null.given.uniform"

model.uniform.prior.plus.traits = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(model.uniform.prior.plus.traits) = rownames(abund)
colnames(model.uniform.prior.plus.traits) = "model.uniform.prior.plus.traits"

model.prior.plus.traits = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(model.prior.plus.traits) = rownames(abund)
colnames(model.prior.plus.traits) = "model.prior.plus.traits"

lambda.values = matrix(0,nrow = length(rownames(abund)),ncol=2+length(rownames(traits)))
rownames(lambda.values) = rownames(abund)
colnames(lambda.values) = c("intercept", rownames(traits), "prior")

delta.R.trait = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(delta.R.trait) = rownames(abund)
colnames(delta.R.trait) = "delta.R.trait"

delta.R.neutral = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(delta.R.neutral) = rownames(abund)
colnames(delta.R.neutral) = "delta.R.neutral"

information.unique.to.local.trait.constraints = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(information.unique.to.local.trait.constraints) = rownames(abund)
colnames(information.unique.to.local.trait.constraints) = "information.unique.to.local.trait.constraints"

information.unique.to.neutral.prior = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(information.unique.to.neutral.prior) = rownames(abund)
colnames(information.unique.to.neutral.prior) = "information.unique.to.neutral.prior"

joint.information = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(joint.information) = rownames(abund)
colnames(joint.information) = "joint.information"

biologically.unexplained.information = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(biologically.unexplained.information) = rownames(abund)
colnames(biologically.unexplained.information) = "biologically.unexplained.information"

pure.trait.effect = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(pure.trait.effect) = rownames(abund)
colnames(pure.trait.effect) = "pure.trait.effect"

pure.metacommunity.effect = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(pure.metacommunity.effect) = rownames(abund)
colnames(pure.metacommunity.effect) = "pure.metacommunity.effect"

joint.trait.metacommunity.effect = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(joint.trait.metacommunity.effect) = rownames(abund)
colnames(joint.trait.metacommunity.effect) = "joint.trait.metacommunity.effect"

unexplained.effect = matrix(0,nrow = length(rownames(abund)),ncol=1)
rownames(unexplained.effect) = rownames(abund)
colnames(unexplained.effect) = "unexplained.effect"

predicted.ra.uniform.plus.traits = matrix(0,nrow = length(rownames(abund)),ncol=length(colnames(abund)))
rownames(predicted.ra.uniform.plus.traits) = rownames(abund)
colnames(predicted.ra.uniform.plus.traits) = colnames(abund)

predicted.ra.prior.plus.traits = matrix(0,nrow = length(rownames(abund)),ncol=length(colnames(abund)))
rownames(predicted.ra.prior.plus.traits) = rownames(abund)
colnames(predicted.ra.prior.plus.traits) = colnames(abund)

#Calculate correlation between local abundances and metacommunity or trait values
cor.local.meta<-cor(as.vector(ra),as.vector(q))
cor.meta.x = cor.local.meta

for(j in 1:length(rownames(traits))){
trait.cor<-cor(as.vector(traits[j,]),as.vector(ra))
cor.trait.x[j] = trait.cor}

#Run first model fit using CWM and only uniform prior
fit1<-maxent2(constr=CWM,states=traits)
temp1<-maxent.test2(model=fit1,obs=ra,nperm=Nperm,quick=FALSE, plot=plot.set)
if(plot.set == T)text(0.001,0.001,"uniform prior")
predicted.ra.uniform.plus.traits = fit1$prob

#An ESTIMATE (using N permutations) of model bias (model.mean.null.given.uniform)
model.mean.null.given.uniform<-temp1$mean.KLR2.null
model.bias<-model.mean.null.given.uniform

#Set model fit using the CWM and only the uniform prior
model.uniform.prior.plus.traits<-temp1$KLR2.prior.plus.traits

#If model using the maximally uninformative prior plus traits is < the model bias, then correct
if(model.uniform.prior.plus.traits<model.bias)model.uniform.prior.plus.traits<-model.bias

#Calculate model using the CWM plus neutral prior
fit2<-maxent2(constr=CWM,states=traits,prior=q,lambda=TRUE)
fitted.lambdas<-fit2$lambda
temp2<-maxent.test2(model=fit2,obs=ra,nperm=Nperm,quick=FALSE, plot=plot.set)
if(plot.set == T)text(0.001,0.001,"neutral prior")
model.mean.null.given.prior<-temp2$mean.KLR2.null
predicted.ra.prior.plus.traits = fit2$prob

#If this null model, given the neutral prior, is < the mean given the maximally uniformative prior, then correct
if(model.mean.null.given.prior<model.bias)model.mean.null.given.prior<-model.bias

#Now fit using the CWM and the neutral prior:
model.prior.plus.traits<-temp2$KLR2.prior.plus.traits

# if this model, given the neutral prior and the traits, is less than the fit of the model using the maximally uniformative prior and the traits, then correct:
if(model.prior.plus.traits<model.uniform.prior.plus.traits)model.prior.plus.traits<-model.uniform.prior.plus.traits

#Now generate the necessary numbers for the output
just.neutral1<-model.prior.plus.traits-model.uniform.prior.plus.traits
just.neutral2<-model.mean.null.given.prior-model.mean.null.given.uniform
just.traits1<-model.uniform.prior.plus.traits-model.mean.null.given.uniform
just.traits2<-model.prior.plus.traits-model.mean.null.given.prior

information.unique.to.local.trait.constraints<-just.traits2
information.unique.to.neutral.prior<-just.neutral1

delta.R.trait<-just.traits1
delta.R.neutral<-just.neutral2
lambda.values[k,]<-fitted.lambdas

T1<-just.traits1-just.traits2
T2<-just.neutral2-just.neutral1

if(abs(T1-T2)>1e-6)return(
 list(message="problemtest",T1=T1,T2=T2,traits1=just.traits1,traits2=just.traits2,neutral1=just.neutral1,neutral2=just.neutral2, diff.traits=abs(just.traits2-just.traits1),diff.neutral=abs(just.neutral1-just.neutral2)))

joint.information<-T1

biologically.unexplained.information<-(1-model.prior.plus.traits)

pure.trait.effect= information.unique.to.local.trait.constraints/(1-model.mean.null.given.uniform) 

pure.metacommunity.effect= information.unique.to.neutral.prior/(1-model.mean.null.given.uniform) 

joint.trait.metacommunity.effect= joint.information/(1-model.mean.null.given.uniform) 

unexplained.effect= biologically.unexplained.information/(1-model.mean.null.given.uniform) 

#Store information on distances
temp_effects = cbind(pure.trait.effect,pure.metacommunity.effect,joint.trait.metacommunity.effect, unexplained.effect)
temp_effects = data.frame(temp_effects)
temp_effects$maxdistance = categories[i]
temp_effects$Plotcode = rownames(m)[k]
temp_effects = temp_effects[,c(6,5,1,2,3,4)]

if (i == 1 & k == 1){ effects_distance <<- temp_effects	} else {
effects_distance <<- rbind(effects_distance,temp_effects)	}

write.table(effects_distance, conn , append = TRUE , col.names = FALSE, row.names = FALSE, sep=",")

if (pure.metacommunity.effect < cutoff){
            stop = TRUE # Fire the flag, and break the inner loop
            break
        }
}  #stop category iteration

if (stop){next}

	} #stop rowname iteration
return(output)
close(conn)
	} #stop foreach 
stopCluster(cl)
		}#stop function





##########################
#### MAXENT2 FUNCTION ####
##########################

maxent2 = function (constr, states, prior, tol = 1e-07, lambda = FALSE)
{
# maxent2 version october 21 2011
 if (is.vector(constr)) {
 means.names <- names(constr)
 constr <- matrix(constr, 1, length(constr))
 dimnames(constr) <- list("set1", means.names)
 }
 if (is.data.frame(constr))
 constr <- as.matrix(constr)
 if (!is.numeric(constr))
 stop("constr must only contain numeric values\n")
 if (!is.numeric(tol))
 stop("tol must be numeric\n")
 if (!is.logical(lambda))
 stop("lambda must be logical\n")
 if (is.vector(states)) {
 s.names <- names(states)
 states <- matrix(states, nrow = 1, ncol = length(states))
 dimnames(states) <- list("constraint", s.names)
 }
 if (is.data.frame(states))
 states <- as.matrix(states)
# new test to exclude negative or zero traits
 if(sum(states<=0)>0)stop("trait values must not be negative or zero")
 if (!is.numeric(states)) 
 stop("states must only contain numeric values\n")
 if (dim(states)[2] == 1 && dim(states)[1] > 1)
 states <- t(states)
 s.names <- dimnames(states)[[2]]
 c.names <- dimnames(states)[[1]]
 set.names <- dimnames(constr)[[1]]
 n.states <- dim(states)[2]
 n.traits <- dim(states)[1]
 n.sets <- dim(constr)[1]
 if (n.traits != dim(constr)[2])
 stop("number of constraints in constr should be equal to number of constraints
in states\n")
 if (missing(prior)) {
 prior <- matrix(1/n.states, n.sets, n.states)
 dimnames(prior) <- list(set.names, s.names)
 }
 if (is.vector(prior)) {
 if (length(prior) != n.states) {
 stop("number of states in prior should be equal to number in states\n")
 }
 if (n.sets == 1) {
 prior <- matrix(prior, 1, length(prior))
 dimnames(prior) <- list("set1", s.names)
 }
 else {
 prior <- matrix(rep(prior, n.sets), n.sets, length(prior),
 byrow = T)
 dimnames(prior) <- list(set.names, s.names)
 }
 }
 if (is.data.frame(prior))
 prior <- as.matrix(prior)
 if (!is.numeric(prior))
 stop("prior must only contain numeric values\n")
 if (dim(prior)[2] == 1 && dim(prior)[1] > 1)
 prior <- t(prior)
 if (dim(prior)[2] != n.states)
 stop("number of states in prior should be equal to number in states\n")
 if (dim(prior)[1] > 1 && dim(prior)[1] != n.sets)
 stop("number of rows in prior should be 1 or equal to number of rows in
constr\n")
 if (dim(prior)[1] == 1)
 prior <- matrix(rep(prior[1, ], n.sets), n.sets, ncol(prior),
 byrow = T)
 if (any(prior < 0 | prior > 1))
 stop("prior must contain probabilities between 0 and 1\n")
 prior <- t(apply(prior, 1, function(x) x/sum(x)))
 dimnames(prior) <- list(set.names, s.names)
 if (any(is.na(constr)) || any(is.na(states)) || any(is.na(prior)))
 stop("no NA's allowed\n")
 allprobs <- matrix(NA, n.sets, n.states)
 dimnames(allprobs) <- list(set.names, s.names)
 moments <- matrix(NA, n.sets, n.traits)
 dimnames(moments) <- list(set.names, c.names)
 entropy <- rep(NA, n.sets)
 names(entropy) <- set.names
 iter <- rep(NA, n.sets)
 names(iter) <- set.names
 if (lambda) {
 lambdas <- matrix(NA, n.sets, n.traits + 2)
 dimnames(lambdas) <- list(set.names, c("intercept", c.names,"prior"))
 }
 for (i in 1:n.sets) {
 itscale <- .Fortran("itscale5", as.double(t(states)),
 as.integer(n.states), as.integer(n.traits), as.double(constr[i,
 ]), as.double(prior[i, ]), prob = double(n.states),
 entropy = double(1), niter = integer(1), as.double(tol),
 moments = double(n.traits), PACKAGE = "FD")
 allprobs[i, ] <- itscale$prob
 moments[i, ] <- itscale$moments
 entropy[i] <- itscale$entropy
 iter[i] <- itscale$niter
 if (lambda)
 lambdas[i, ] <- coef(lm(log(itscale$prob) ~ t(states)+log(prior[i,])))
 }
 res <- list()
 if (n.sets == 1) {
 res$prob <- allprobs[1, ]
 res$moments <- moments[1, ]
 names(entropy) <- NULL
 res$entropy <- entropy
 names(iter) <- NULL
 res$iter <- iter
 if (lambda)
 res$lambda <- lambdas[1, ]
 res$constr <- constr[1, ]
 res$states <- states
 res$prior <- prior[1, ]
 }
 else {
 res$prob <- allprobs
 res$moments <- moments
 res$entropy <- entropy
 res$iter <- iter
 if (lambda)
 res$lambda <- lambdas
 res$constr <- constr
 res$states <- states
 res$prior <- prior
 }
 return(res)
}




###############################
#### MAXENT2.TEST FUNCTION ####
###############################

maxent.test2 = function (model, obs, sub.c, nperm = 999, quick = TRUE, alpha = 0.05, plot = TRUE, minperms = 20)
{
# maxent.test2 version 21 october 2011
 if (is.vector(obs)) {
 s.names <- names(obs)
 obs <- matrix(obs, 1, length(obs))
 dimnames(obs) <- list("set1", s.names)
 }
 if (is.data.frame(obs))
 obs <- as.matrix(obs)
 obs.names <- dimnames(obs)
 if (!is.numeric(obs))
 stop("obs must only contain numeric values\n")
 if (dim(obs)[2] == 1 && dim(obs)[1] > 1)
 obs <- t(obs)
 if (any(is.na(obs)))
 stop("no NA's allowed\n")
 if (any(obs < 0 | obs > 1))
 stop("obs must contain probabilities between 0 and 1\n")
 obs <- t(apply(obs, 1, function(x) x/sum(x)))
 dimnames(obs) <- obs.names
 states <- model$states
# new test to exclude negative or zero trait values...
 if(sum(states<=0)>0)stop("trait values must not be negative or zero")
 s.names <- dimnames(states)[[2]]
 constr <- model$constr
 prob <- model$prob
 prior <- model$prior
# n.states is the # species in the species pool
 n.states <- dim(states)[2]
# 17 nov change
# This is the generic KLR2 function to calculate fit against
# the uniform distribution
 KLR2.fit <- function(oi, pi, n.states){
 sel<-oi>0
 qi.uniform<-rep(1/n.states,length(oi))
 1-sum(oi[sel]*log(oi[sel]/pi[sel]))/
 sum(oi[sel]*log(oi[sel]/qi.uniform[sel]))
 }
# 17 nov change
# observed KLR2...
 sel<-obs>0
 oi<-as.double(as.vector(as.matrix(obs[sel])))
 pi<-as.double(as.vector(as.matrix(prob[sel])))
 KLR2.prior.plus.traits<-KLR2.fit(oi,pi,n.states)
#......
 if (is.vector(constr)) {
 means.names <- names(constr)
 constr <- matrix(constr, 1, length(constr))
 dimnames(constr) <- list("set1", means.names)
 prior <- matrix(prior, 1, length(prior))
 dimnames(prior) <- list("set1", s.names)
 prob <- matrix(prob, 1, length(prob))
 dimnames(prob) <- list("set1", s.names)
 }
# n.sets is # of local communities
 n.sets <- dim(constr)[1]
 if (n.sets != dim(obs)[1])
 stop("number of rows in obs and constr should be equal\n")
 stat <- function(o, p, q) sum(o * log(p/q))
 values <- rep(NA, nperm)
 KLR2.null <- rep(NA, nperm)
 oi <- as.double(as.vector(as.matrix(obs)))
 qi<-as.double(as.vector(as.matrix(prior)))
 sel.oi <- (oi > 0)
# next section used if testing all traits together, not a subset
 if (missing(sub.c)) {
 obs.stat <- stat(obs, prob, prior)
 count <- 0
 outside.ci <- FALSE
 while ((!outside.ci && count < nperm) | count < minperms) {
 count <- count + 1
 prob.temp <- matrix(NA, n.sets, n.states)
 for (j in 1:n.sets) {
 shuffled <- sample(1:n.states, n.states)
 states.perm <- states[, shuffled, drop = F]
 colnames(states.perm) <- s.names
 constr.perm <- functcomp(t(states.perm), obs[j,
 , drop = F])
# maxent2 replaces the old maxent function
# prob.temp[j, ] <- maxent(constr.perm, states.perm,
 prob.temp[j, ] <- maxent2(constr.perm, states.perm,
 prior[j, ])$prob
 }
 pi <- as.double(as.vector(as.matrix(prob.temp)))
# KLR2.null still includes information from the prior (if not uniform)
# 17 nov
KLR2.null[count]<-KLR2.fit(oi,pi,n.states)
# KLR2.null[count] <- 1 - sum(oi[sel.oi] * log(oi[sel.oi]/pi[sel.oi]))/sum(oi[sel.oi] *
# log(oi[sel.oi]/qi.[sel.oi]))
 values[count] <- stat(obs, prob.temp, prior)
 val.temp <- values[1:count]
 p.temp <- (length(val.temp[val.temp >= obs.stat]) +
 1)/(length(val.temp) + 1)
 if (quick && p.temp < 1) {
 ci.hi <- p.temp + 1.96 * sqrt((p.temp * (1 -
 p.temp))/count)
 ci.lo <- p.temp - 1.96 * sqrt((p.temp * (1 -
 p.temp))/count)
 outside.ci <- ci.hi <= alpha || ci.lo >= alpha
 if (outside.ci && count > minperms)
 nperm = count
 }
 }
 }
# next section is if only a subset of constraints are being tested.
 else {
 if (length(sub.c) >= n.states)
 stop("sub.c contains as many or more elements than there are states\n")
 if (is.character(sub.c) && !all(sub.c %in% dimnames(states)[[1]]))
 stop("sub.c does not match the names of the state attributes\n")
 if (is.character(sub.c) && !all(sub.c %in% dimnames(constr)[[2]]))
 stop("sub.c does not match the constraint names\n")
 if (is.character(sub.c))
 sub.c <- which(dimnames(states)[[1]] %in% sub.c)
 if (!is.vector(sub.c))
 stop("sub.c must be a vector\n")
 prob.a <- maxent(constr[, -sub.c, drop = F], states[-sub.c,
 , drop = F], prior)$prob
 obs.stat <- stat(obs, prob, prob.a)
 count <- 0
 outside.ci <- FALSE
 while ((!outside.ci && count < nperm) | count < minperms) {
 count <- count + 1
 prob.temp <- matrix(NA, n.sets, n.states)
 for (j in 1:n.sets) {
 shuffled <- sample(1:n.states, n.states)
 states.perm <- states
 states.perm[sub.c, ] <- states[sub.c, shuffled,
 drop = F]
 colnames(states.perm) <- s.names
 constr.perm <- functcomp(t(states.perm), obs[j,
 , drop = F])
# maxent2 replaces the old maxent
# prob.temp[j, ] <- maxent(constr.perm, states.perm,
 prob.temp[j, ] <- maxent2(constr.perm, states.perm,
 prior[j, ])$prob
 }
 pi <- as.double(as.vector(as.matrix(prob.temp)))
# KLR2.null still includes information from the prior (if not uniform)
# 17 nov
KLR2.null[count]<-KLR2.fit(oi,pi,n.states)
# KLR2.null[count] <- 1 - sum(oi[sel.oi] * log(oi[sel.oi]/pi[sel.oi]))/sum(oi[sel.oi] * log(oi[sel.oi]/qi[sel.oi]))
 values[count] <- stat(obs, prob.temp, prob.a)
 val.temp <- values[1:count]
 p.temp <- (length(val.temp[val.temp >= obs.stat]) +
 1)/(length(val.temp) + 1)
 if (quick && p.temp < 1) {
 ci.hi <- p.temp + 1.96 * sqrt((p.temp * (1 -
 p.temp))/count)
 ci.lo <- p.temp - 1.96 * sqrt((p.temp * (1 -
 p.temp))/count)
 outside.ci <- ci.hi <= alpha || ci.lo >= alpha
 if (outside.ci && count > minperms)
 nperm = count
 }
 }
 }
 values <- values[!is.na(values)]
 p.perm <- (length(values[values >= obs.stat]) + 1)/(length(values) +
 1)
 p.perm.hi <- p.perm + 1.96 * sqrt((p.perm * (1 - p.perm))/nperm)
 p.perm.lo <- p.perm - 1.96 * sqrt((p.perm * (1 - p.perm))/nperm)
 opqfit <- function(o, p, q) 1 - (sum((o - p)^2)/sum((o -
 q)^2))
 r2.op <- cor(as.double(prob), as.double(obs))^2
 if (length(unique(as.double(prior))) != 1)
 r2.oq <- cor(as.double(prior), as.double(obs))^2
 else r2.oq = 0
 fit <- opqfit(obs, prob, prior)
 if (!missing(sub.c)) {
 r2.opa <- cor(as.double(prob.a), as.double(obs))^2
 fit.a <- opqfit(obs, prob.a, prior)
 }
 mean.KLR2.null<-mean(KLR2.null[!is.na(KLR2.null)])
 if (plot) {
 if (missing(sub.c)) {
 par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
 fit.text <- bquote(fit[bold(o * "," * p) * "|" *
 bold(q)] == .(round(fit, 3)))
 r2.op.text <- bquote(italic(KLR)^2 == .(round(KLR2.prior.plus.traits,
 3)))
 r2.oq.text <- bquote(italic(KLR)^2 == .(round(mean.KLR2.null,
 3)))
 plot(as.double(obs), as.double(prob), xlim = c(0,
 1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities",
 main = expression(bold(p)(bold(C) * "," * ~bold(q))))
 abline(0, 1, col = "grey25")
 text(0.1, 0.9, fit.text, cex = 1.2, pos = 4)
 text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4)
 lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed",
 col = "grey50")
 lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed",
 col = "grey50")
 mtext("arithmetic scale", line = 0)
 plot(as.double(obs) + 1e-04, as.double(prob) + 1e-04,
 xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed
probabilities",
 ylab = "predicted probabilities", log = "xy",
 main = expression(bold(p)(bold(C) * "," * ~bold(q))))
 abline(0, 1, col = "grey25")
 lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04,
 1 + 1e-04), lty = "dashed", col = "grey50")
 lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 +
 1e-04), lty = "dashed", col = "grey50")
 mtext("log10 scale, + 1e-4", line = 0)
 plot(as.double(obs), as.double(prior), xlim = c(0,
 1), ylim = c(0, 1), xlab = "", ylab = "prior",
 main = expression(bold(q)))
 abline(0, 1, col = "grey25")
 text(0.1, 0.9, r2.oq.text, cex = 1.2, pos = 4)
 lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed",
 col = "grey50")
 lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed",
 col = "grey50")
 mtext("arithmetic scale", line = 0)
 plot(as.double(obs) + 1e-04, as.double(prior) + 1e-04,
 xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed
probabilities",
 ylab = "prior", log = "xy", main = expression(bold(q)))
 abline(0, 1, col = "grey25")
 lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04,
 1 + 1e-04), lty = "dashed", col = "grey50")
 lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 +
 1e-04), lty = "dashed", col = "grey50")
 mtext("log10 scale, + 1e-4", line = 0)
 mtext(expression(H[0] %->% bold(p)(bold(C) * "," *
 ~bold(q)) == bold(q)), cex = 1.2, outer = T,
 las = 1, line = 1)
# here
 iter.text<-bquote(italic(Nperms) == .(round(count,0)))
 p.text <- bquote(italic(P) == .(round(p.perm, 3)))
 mtext(p.text, cex = 1.2, outer = T, las = 1, line = -0.5)
 mtext(iter.text, cex=0.75, outer=T, las=1, line=-1.5)
 }
 else {
 par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
 fit.text <- bquote(fit[bold(o * "," * p) * "|" *
 bold(q)] == .(round(fit, 3)))
 fit.a.text <- bquote(fit[bold(o * "," * p) * "|" *
 bold(q)] == .(round(fit.a, 3)))
 r2.op.text <- bquote(italic(KLR)^2 == .(round(KLR2.prior.plus.traits,
 3)))
 r2.opa.text <- bquote(italic(KLR)^2 == .(round(KLR2.null,
 3)))
 plot(as.double(obs), as.double(prob), xlim = c(0,
 1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities",
 main = expression(bold(p)(bold(A) ~ union(bold(B) *
 "," * ~bold(q)))))
 abline(0, 1, col = "grey25")
 text(0.1, 0.9, fit.text, cex = 1.2, pos = 4)
 text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4)
 lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed",
 col = "grey50")
 lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed",
 col = "grey50")
 mtext("arithmetic scale", line = 0)
 plot(as.double(obs) + 1e-04, as.double(prob) + 1e-04,
 xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed
probabilities",
 ylab = "predicted probabilities", log = "xy",
 main = expression(bold(p)(bold(A) ~ union(bold(B) *
 "," * ~bold(q)))))
 abline(0, 1, col = "grey25")
 lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04,
 1 + 1e-04), lty = "dashed", col = "grey50")
 lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 +
 1e-04), lty = "dashed", col = "grey50")
 mtext("log10 scale, + 1e-4", line = 0)
 plot(as.double(obs), as.double(prob.a), xlim = c(0,
 1), ylim = c(0, 1), xlab = "", ylab = "", main =
expression(bold(p)(bold(A) *
 "," * ~bold(q))))
 abline(0, 1, col = "grey25")
 text(0.1, 0.9, fit.a.text, cex = 1.2, pos = 4)
 text(0.1, 0.75, r2.opa.text, cex = 1.2, pos = 4)
 lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed",
 col = "grey50")
 lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed",
 col = "grey50")
 mtext("arithmetic scale", line = 0)
 plot(as.double(obs) + 1e-04, as.double(prob.a) +
 1e-04, xlim = c(1e-04, 1), ylim = c(1e-04, 1),
 xlab = "observed probabilities", ylab = "", log = "xy",
 main = expression(bold(p)(bold(A) * "," * ~bold(q))))
 abline(0, 1, col = "grey25")
 lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04,
 1 + 1e-04), lty = "dashed", col = "grey50")
 lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 +
 1e-04), lty = "dashed", col = "grey50")
 mtext("log10 scale, + 1e-4", line = 0)
 mtext(expression(H[0] %->% bold(p)(bold(A) ~ union(bold(B)) *
 "," * ~bold(q)) == bold(p)(bold(A) * "," * ~bold(q))), 
 cex = 1.2, outer = T, las = 1, line = 1)
 iter.text<-bquote(italic(Nperms) == .(round(count,0)))
 p.text <- bquote(italic(P) == .(round(p.perm, 3)))
 mtext(p.text, cex = 1.2, outer = T, las = 1, line = -0.5)
 mtext(iter.text, cex=0.75, outer=T, las=1, line=-1.5)
 }
 }
 res <- list()
 res$fit <- fit
 if (!missing(sub.c)) {
 res$fit.a <- fit.a
 res$r2.a <- r2.opa
 }
# res$r2 <- r2.op
# if (missing(sub.c))
# res$r2.q <- r2.oq
 res$obs.stat <- obs.stat
 res$nperm <- nperm
 res$pval <- p.perm
 res$ci.pval <- c(p.perm.lo, p.perm.hi)
 res$KLR2.null<-KLR2.null[!is.na(KLR2.null)]
 res$mean.KLR2.null <- mean(KLR2.null[!is.na(KLR2.null)])
 res$KLR2.prior.plus.traits<-KLR2.prior.plus.traits
 return(res)
}
