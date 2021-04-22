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

#!!!!!help dont know what trait data file is for our data!!!
#maybe samples.txt?

### Outlier detection incorporated into trait measures. 
traitData= read.csv("samples.csv", row.names=1)
dim(traitData)
head(traitData)
names(traitData)

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

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed")
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

save(datExpr0, datTraits, file="Crep_Samples_Traits_ALL.RData")

################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() use this in base R
allowWGCNAThreads() 
lnames = load(file="Crep_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust these power values to hone in on proper sft value
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=13 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

minModuleSize=90 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

# # dynamicMods
#dynamicMods
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19 
#1164  576  415  325  272  263  259  258  249  244  241  236  232  222  219  212  206  198  192 
#20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38 
#185  184  183  182  180  178  177  176  175  172  170  165  163  158  158  157  157  153  152 
#39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54   55   56   57 
#150  145  144  141  137  137  135  134  133  133  131  130  129  126  123  118  118  116  112 
#58   59   60   61   62   63   64   65   66 
#111  104  100   99   97   96   94   93   91

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 13)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_crep_nomerge.RData")

lnames = load(file = "Network_crep_nomerge.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.6
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.6.RData")

###############Relating modules to traits and finding important genes
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Crep_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.6.RData");
lnames = load(file = "Network_crep_nomerge.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

#moduleColors
#black         blue       coral2     darkgrey  floralwhite       grey60        ivory 
#1338         2928         1167         1088          895          206          564 
#mediumorchid   orangered4       purple    steelblue 
#478         1557          615         1749 

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
quartz()
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

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$pH7.5); #change Lipidrobust to your trait name
names(weight) = "pH7.5"
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
module = "coral2"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for pH7.5",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="coral2"]) #black  blue brown green  grey  pink   red 

c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors)
#moduleColors
#black         blue       coral2     darkgrey  floralwhite       grey60        ivory 
#1338         2928         1167         1088          895          206          564 
#mediumorchid   orangered4       purple    steelblue 
#478         1557          615         1749
head(c.vsd)
write.csv(c.vsd,"rlog_MMcoral2.csv",quote=F)

##############################heatmap of module expression with bar plot of eigengene, no resorting of samples...
#names(dis)
sizeGrWindow(8,7);
which.module="coral2" #pick module of interest
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

##############################heatmap of module expression with bar plot of trait of interest by sample...
#here we just have binary traits, but if you have a continuous trait this code is cool
sizeGrWindow(8,7);
which.module="yellow" #pick module of interest
which.trait="fvfm" #change trait of interest here
datTraits=datTraits[order((datTraits$fvfm),decreasing=T),]#change trait of interest here

trait=datTraits[, paste(which.trait)]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="fvfm",xlab="sample")#change trait of interest here

#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMcoral2.csv", row.names=1)
head(vsd)
gg=read.table("Crep454_iso2gene.tab", sep="\t")
head(gg)
library(pheatmap)

############################################
whichModule="coral2"
top=100

datME=MEs
vsd <- read.csv("Crep_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)

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

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))

###fisher for GO
##########fisher of module vs whole dataset
library(WGCNA)
vsd <- read.csv("Crep_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="coral2" # name your color and execute to the end

length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_fisher.csv",sep=""),quote=F)

modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

######--------------------end--------------------#######
