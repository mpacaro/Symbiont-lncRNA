library(WGCNA)
library(flashClust)


options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("lcpm.csv") #bringing in filtered and normalized count data from edgeR
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#32095, this is the same as before (32095 transcripts/isoforms - about 2/3 (67%) of the number that we started with)

datExpr0 = as.data.frame(t(dat))
View(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes, if false run the script below
#this is true for us, so we don't need to run the code below

#if (!gsg$allOK)
#{if (sum(!gsg$goodGenes)>0)
  #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
  #if (sum(!gsg$goodSamples)>0)
    #printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  #datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#}
#gsg=goodSamplesGenes(datExpr0, verbose = 3)
#gsg$allOK 
#dim(datExpr0) 


### Outlier detection incorporated into trait measures.
setwd("/projectnb/bi594/mpacaro/lncRNA/")
getwd()
traitData= read.csv("samples_traits.csv", row.names=1)
dim(traitData)
head(traitData)
names(traitData) #check to make sure names are the same (traits table and gene expression files)

# Form a data frame analogous to expression data that will hold the clinical traits.
dim(datExpr0)
rownames(datExpr0)
rownames(traitData)=rownames(datExpr0)
traitData$Sample= NULL 
# datTraits=allTraits
datTraits=traitData

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#making adjaceny matrix to show how linked genes are
#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed") #default is unsigned, we are using signed here
# this calculates the whole network connectivity we choose signed because we care about direction of gene expression
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

#looking at overall similarity of sames and which traits they have 
#code above is making this dendrogram


# Remove outlying samples from expression and trait data
# remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# datExpr=datExpr0[!remove.samples,]
# datTraits=datTraits[!remove.samples,]

save(datExpr0, datTraits, file="lncRNA_Samples_Traits_ALL.RData")
#save what we just did as Rdata


################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() use this in base R
allowWGCNAThreads() 
lnames = load(file="lncRNA_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust these power values to hone in on proper sft value
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)
sft
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"), ylim=c(0,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#making these plots allows us to pick our soft threshold - pick the value that is under the red line
#here our R^2 values do not go above the 0.90 line.. we are continuing and picking softpower=15 bc thats where it looks like it levels off
#mary also used 0.90

softPower=15 #smallest value to plateau at ~0.80
adjacency=adjacency(datExpr0, power=softPower,type="signed") 
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

#this part above takes a while

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
#ridiculogram - y axis is dissimilarity so lines that extend far down are genes most closely related, genes at the top are more dissimilar

#chosing module size, will be diff for each proj
#module size = number of genes in each module, they picked 90 for GO enrichment downstream
minModuleSize=90 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)
#should show modules and they should all be bigger than 90
#each module gets number and size 

#dynamicMods
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#6730 6224 5578 5184 1453 1147 1011  935  838  668  539  536  505  395  352 

#now we are assigning colors instead of the numbers above
dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 15)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_lncrna_nomerge.RData")

lnames = load(file = "Network_lncrna_nomerge.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "") #save this figure and showwith modtrait heatmap at same time

MEDissThres= 0.35  #0.6 #start with 0, look at modtrait heatmap
#0.35 gives 5 modules

abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#ordering our colors here 
moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.35.RData")

###############Relating modules to traits and finding important genes
#bring in our traits and bring in 2 datasets 
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "lncRNA_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.

lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.35.RData");
lnames = load(file = "Network_lncrna_nomerge.RData");
lnames

#how wgcna calls our modules
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

#moduleColors for 0.35 threshold
#brown greenyellow     magenta      salmon   turquoise 
#14157         891        2441         505       14101 

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
#quartz()
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#include trait heatmap
#the number in the heatmap can be interpreted as ___% of the variation is explained by this trait

#for example, can look into highly correlated 

#potential modules of interest
#BROWN - high correlation with sensitive.temp27
#TURQUOISE - high correlation with sensitive.temp27

#GREENYELLOW - high correlation with temp27 and temp32
#MAGENTA - high correlation with temp 27 and temp 32

#SALMON - high correlation with sensitive and tolerant 


#i think look at magenta with each of the temps and then brown and turqouise at sensitive. temp 27

#######################this is code from Mary's script to find # of mRNA and lncRNA in each module##################################
# genes found in each module 
gyell <- names(datExpr0)[moduleColors=="greenyellow"] 
turq <- names(datExpr0)[moduleColors=="turquoise"] 
mag <- names(datExpr0)[moduleColors=="magenta"] 
brwn <- names(datExpr0)[moduleColors=="brown"] 
salmon <- names(datExpr0)[moduleColors=="salmon"] 
# lncNRAs
# number of lncRNAs in each module 
length(grep("^MSTRG", gyell)) #118
length(grep("^MSTRG", turq)) #2271
length(grep("^MSTRG", mag)) #300
length(grep("^MSTRG", brwn)) #5149
length(grep("^MSTRG", salmon)) #59

#sum(118+2271+300+5149+59), sanity check, this matches the total length of lncRNAs below
## merge the lncRNAs from each module into one big file   
lncRNA_allmodules <- c(gyell[grep("^MSTRG", gyell)], turq[grep("^MSTRG", turq)], mag[grep("^MSTRG", mag)], brwn[grep("^MSTRG", brwn)], salmon[grep("^MSTRG", salmon)])
length(lncRNA_allmodules) # 7897
## remove any repeated genes 
uniq_lncRNA_allmodules <- unique(lncRNA_allmodules)
length(uniq_lncRNA_allmodules)
# 7897 lncRNAs - no overlapping lncRNAs among modules 

# mRNAs
# number of mRNAs in each module 
length(grep("^mRNA", gyell)) #773
length(grep("^mRNA", turq)) #11830
length(grep("^mRNA", mag)) #2141
length(grep("^mRNA", brwn)) #9008
length(grep("^mRNA", salmon)) #446

## merge the mRNAs from each module into one big file   
mRNA_allmodules <- c(gyell[grep("^mRNA", gyell)], turq[grep("^mRNA", turq)], mag[grep("^mRNA", mag)], brwn[grep("^mRNA", brwn)], salmon[grep("^mRNA", salmon)])
length(mRNA_allmodules) #24198
## remove any repeated genes 
uniq_mRNA_allmodules <- unique(mRNA_allmodules)
length(uniq_mRNA_allmodules)
# 24198 mRNAs - no overlapping mRNAs among modules 
# all of the lncRNAs and mRNAs fell in one of the modules and the results were different from the Pearson correlation results possibly because in the Pearson correlation, we had the threshold cutoffs for the coefficient and p-value which could have lowered the number of lncRNAs while in WGCNA all of the genes are assigned to a module without any cutoff

###############################################

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$temp32); #change Lipidrobust to your trait name

#not sure what the trait name should be 
names(weight) = "temp32"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Gene-trait significance correlation plots
# par(mfrow=c(2,3))
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for temp32",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#this is correlation plot ^^^
#brown module and sensitive.temp27 are highly correlated
#brown module and sensitive.temp32 not that correlated
#could create these figures for each of the modules you find interesting that you want to further explore 


#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="magenta"]) 


#making the figure below look wrong but i think its just bc we have a lot of samples?
c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size,  brown = 14157
table(moduleColors)
#moduleColors
#brown greenyellow     magenta      salmon 
#14157         891        2441         505 
#turquoise 
#14101 
head(c.vsd)
write.csv(c.vsd,"rlog_MMmagenta.csv",quote=F)
#this is all of the genes in the brown module it is creating csv with that subset of data



##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="magenta" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")
#this is a cool plot where you can see that genes in this module are upregulated in the pH7.5 treatment
#this is looking at eigengene expression (overall expression of where genes are going) up = upreg, down=downreg
#~1:18

##############################heatmap of module expression with bar plot of trait of interest by sample...
#here we just have binary traits, but if you have a continuous trait this code is cool
#sizeGrWindow(8,7);
#which.module="yellow" #pick module of interest
##which.trait="fvfm" #change trait of interest here
#datTraits=datTraits[order((datTraits$fvfm),decreasing=T),]#change trait of interest here

#trait=datTraits[, paste(which.trait)]
#genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
#genes=genes[rownames(datTraits),]

#quartz()
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
#par(mar=c(5, 4.2, 0, 0.7))
#barplot(trait, col=which.module, main="", cex.main=2,
        #ylab="fvfm",xlab="sample")#change trait of interest here

#getting error here !!!!!!!!!!!!!!!!!!!!!!!!!!!

#this is KME data
#Gene relationship to trait and important modules: Gene Significance and Module membership


allkME =as.data.frame(signedKME(datExpr0, MEs)) #changed t(dat) to datExpr0 because problems with extra X row 
head(allkME)
vsd=read.csv(file="rlog_MMmagenta.csv", row.names=1)
head(vsd)
gg=read.table("transcript2geneDescription.tab", sep="\t") 
head(gg)
library(pheatmap)

############################################
whichModule="magenta"
top=100 #looking at the top 100 genes 

datME=MEs
vsd <- read.csv("lcpm.csv", row.names=1) 
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
summary(hubs)

#labeling row names with g names from gg
gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames
length(hubs)
#48 that are significantly differentially expressed??


contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
sizeGrWindow(12,9)
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))
#makes heatmap with annotations of genes if they have them
#genes that assign to color module the best are shown here 
#~1:21
#this window is really bad and hard to see what this heatmap says??????????????????

###fisher for GO - 1 if in module, 0 if not 
#this is binary way
#can do this way or KME way (more for continous data?)
##########fisher of module vs whole dataset
library(WGCNA)
vsd <- read.csv("lcpm.csv", row.names=1) 
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
View(data)
allkME =as.data.frame(signedKME(data, MEs))


whichModule="magenta" # name your color and execute to the end

#now it is just asking which genes are in the module and going to give a 1 or 0
length(moduleColors) 
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1]) #should be same number of color module chosen above, sanity check
head(inModule)
View(inModule)
write.csv(inModule,file=paste(whichModule,"_fisher.csv",sep=""),quote=F) #this will automatically add the color name you are in to csv
#will be file with samples and column of 0 or 1


#i think this is for continuous data and not needed for us
#KME way is how good gene fits to module?
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

######--------------------end--------------------#######
