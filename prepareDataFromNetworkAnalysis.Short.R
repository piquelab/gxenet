#######################################################################
## Extract RPKM and build network 
## 1) Count to FPKMs
## 2) Filter out genes with <10FPKM for >90% of columns 
## 3) Remove GC content effect 
## 4) Quantile Normalization
## 5) Adjust for individual effect across samples (mean)
## 6) One step Automated Network construction (soft threshold power 12)
#########################################################################


##Load Modules
#sink(file="stdout.txt",type=c("output","message"))
library(WGCNA)

##library(biomaRt)
library(gplots)
library(cluster)
library(parallel)
library(gtools)
LPG <- Sys.getenv("LPG")
options(stringsAsFactors = FALSE);

load("~/piquelab/Allison/SystemPaper/network/CytoscapeFiles/NetworkEverything.Rd")

treatment="T6C1";
celltype="HUVEC";
threshold=0.05;
OfInterestMods=72;


## treatmentKEY
## allcv2  covariate table, may not need.
## comb combinations parsed allcv2

names(comb)<- c("tr","ct","value")


cytlist <- lapply(rownames(comb),function(x){
  lapply(1:87,function(m){
  tr <- as.character(comb[x,"tr"])
  ct <- as.character(comb[x,"ct"])
  cat("Processing:",x,tr,ct,m,"\n")
  cyto <- Vizmodule(m,tr,threshold,ct)
  })
})

save(comb,allcv2,treatmentKEY,cytlist,file="net.Rd")



Vizmodule <- function(OfInterestMods, treatment, threshold, celltype)
{
	inModule <- is.finite(match(net$colors,OfInterestMods))
	modTranscripts <-transcripts[inModule]
	modGenes <- ExprTbl[inModule,]$g.id
	modTOM <- asM.TOM[inModule, inModule]
	dimnames(modTOM) <- list(modTranscripts, modTranscripts)
        DEdata <- read.table(allDEfiles[celltype][[1]][treatment],header=T,as.is=T,sep=" ",row.names="t.id")
        DEdata <- DEdata[rownames(DEdata) %in% rownames(get(paste0(celltype,"datatrans"))),]
        nodeAttr <- DEdata[DEdata$ensg %in% rownames(ExprTbl),]
	#DExData <- as.data.frame(DEdfs[treatment])
        #colnames(DExData) <- sapply(strsplit(colnames(DExData),".",fixed=T),'[',2)
        #Now we add data of interest to the Node Attribute data.frame, Including DEGs pvalues, Treatment specific t-Values,Connectivity etc. 
	nodeAttr <- DEdata[match(transcripts[inModule],DEdata$ensg),]
        nodeAttr$DegAdjPval <- (nodeAttr$padj <= 0.01) * nodeAttr$logFC
        IDX <- as.integer(moduleLabels) %in% as.integer(OfInterestMods)
	nodeAttr$ModConnec <- ExprTbl$kWithin[IDX] #Degree within Module
	nodeAttr$NetConnec <- ExprTbl$kTotal[IDX]  #Degree within Network
        nodeAttr$TreatmentColor <- ExprTbl$TreatmentColor[IDX] #Treatment Colors
        nodeAttr$maxT <- ExprTbl$maxT[IDX] #Maximum T value
        nodeAttr$mods <- ExprTbl$mods[IDX]
	#nodeAttr$CorConnec <- abs(dataKME[IDX ,as.numeric(OfInterestMods)])  #Correlation to module Eigengene
        nodeAttr$test <- ExprTbl$test[IDX]
        nodeAttr$ZscoreCols <- ExprTbl$ZscoreCols[IDX]
        nodeAttr$NodeSig <- -log2(nodeAttr$padj)
        nodeAttr$Specificity <- ExprTbl$Specificity[IDX]
        ### Adding color to node
        breaks <- c(-10,((0:201/101)-1)*1.0,10) 
        bin <- cut(nodeAttr$logFC,breaks)
        lev <- levels(bin)
        lev.col <- bluered(length(lev))
        nodeAttr$logFC.color <- lev.col[bin]
        nodeAttr$logFC.color[is.na(nodeAttr$logFC.color)] = "#FFFFFF"
        ##
        cyt <- nodeAttr
        cyt <- exportCytoscape(modTOM,
		edgeFile = paste("./cyt/CytoscapeInput-edges-", paste0(OfInterestMods[1],"-", treatment,"-",celltype), ".txt", sep=""),
		nodeFile = paste("./cyt/CytoscapeInput-nodes-", paste0(OfInterestMods[1],"-", treatment,"-",celltype), ".txt", sep=""),
		weighted = TRUE,
		threshold = threshold,
		nodeNames = modGenes,
		altNodeNames= modTranscripts,
		nodeAttr = nodeAttr);
	return(cyt)
}


























source("/wsu/home/fs/fs92/fs9227/piquelab/charvey/GxE_041015/network_analysis/scripts/heatmap.multi.R")
##Set working directory
#setwd("/wsu/home/fs/fs92/fs9227/network_analysis/")

celline <- "HUVEC"
cellines <- c("HUVEC","LCL","Mel","SMC","PBMC")
plateKey <- data.frame(c("DP1","DP7"))
colnames(plateKey) <- "LCL"
plateKey$SMC <- c("DP10","DP9")
plateKey$Mel <- c("DP11","DP12")
plateKey$PBMC <- c("DP5","DP6")
plateKey$HUVEC <- c("DP4","DP8")
platePrefix <- "DP4"           #Default value
platePrefixs <- c("DP4","DP8") #Default Values



#Load Treatmet Color data.frame
ColKey <- read.csv("~/Network_Analysis/ColorKey.csv", header=T)

##Define Functions
#Parallel sapply
cores <- as.integer(Sys.getenv("NCPUS"))
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefixs <- c(cargs[1],cargs[2])
  celline <- cargs[3]
  if(length(cargs)>=4)
  cores <- as.numeric(cargs[4])


if(cores<1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}


#Allow multiple cores on Network Build
enableWGCNAThreads(cores)

#Export cytoscape network edge and node files while allowing nodeAttr to be a data.frame

exportCytoscape <- function (adjMat, edgeFile = NULL, nodeFile = NULL, weighted = TRUE,
    threshold = 0.5, nodeNames = NULL, altNodeNames = NULL, nodeAttr = NULL,
    includeColNames = TRUE)
{
    adjMat = as.matrix(adjMat)
    adjMat[is.na(adjMat)] = 0
    nRow = nrow(adjMat)
    checkAdjMat(adjMat, min = -1, max = 1)
    if (is.null(nodeNames))
        nodeNames = dimnames(adjMat)[[1]]
    if (is.null(nodeNames))
        stop("Cannot determine node names: nodeNames is NULL and adjMat has no dimnames.")
    rowMat = matrix(c(1:nRow), nRow, nRow, byrow = TRUE)
    colMat = matrix(c(1:nRow), nRow, nRow)
    adjDst = as.dist(adjMat)
    dstRows = as.dist(rowMat)
    dstCols = as.dist(colMat)
    edges = abs(adjDst) > threshold
    nEdges = sum(edges)
    edgeData = data.frame(fromNode = nodeNames[dstRows[edges]],
        toNode = nodeNames[dstCols[edges]], weight = if (weighted)
            adjDst[edges]
        else rep(1, nEdges), direction = rep("undirected", nEdges),
        fromAltName = if (is.null(altNodeNames))
            rep("NA", nEdges)
        else altNodeNames[dstRows[edges]], toAltName = if (is.null(altNodeNames))
            rep("NA", nEdges)
        else altNodeNames[dstCols[edges]])
    nodesPoresent = rep(FALSE, ncol(adjMat))
    #print("Here2")
    nodesPresent[dstRows[edges]] = TRUE
    nodesPresent[dstCols[edges]] = TRUE
    nNodes = sum(nodesPresent)
    #print("Here3")
    nodeData = data.frame(nodeName = nodeNames[nodesPresent],
        altName = if (is.null(altNodeNames))
            rep("NA", nNodes)
        else altNodeNames[nodesPresent], nodeAttribute = if (is.null(nodeAttr))
            rep("NA", nNodes)
        else nodeAttr[nodesPresent,])
    if (!is.null(edgeFile))
        write.table(edgeData, file = edgeFile, quote = FALSE,
            row.names = FALSE, col.names = includeColNames)
    if (!is.null(nodeFile))
        print(head(nodeData))
        write.table(nodeData, file = nodeFile, quote = FALSE,
            row.names = FALSE, col.names = includeColNames)
    list(edgeData = edgeData, nodeData = nodeData)
}


#Generate cytoscape ouput for a specific module and the treatment it responds to and the TOM threshold


plot.multi.dens <- function(s,xlabel,ylabel,cols)
  {
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
      junk.x = c(junk.x, density(s[[i]])$x)
      junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = "",xlab=xlabel,ylab=ylabel)
    for(i in 1:length(s)) {
      lines(density(s[[i]]), xlim = xr, ylim = yr, col = cols[i])
    }
  }



##Load FPKM data from 2 plates


load(file=paste0("/wsu/home/groups/piquelab/charvey/GxE_041015/network_analysis/Full/fpkms/",platePrefixs[1],".counts.fpkm.gc.qn.adj.Rd"))
P1cv2 <- cv2
P1fpkms <- fpkms
rm(list=c("fpkms","cv2")) #Cleanup

#allfpkms <- merge(P1fpkms,allfpkms, by=0, all=FALSE)
#rownames(allfpkms) <- allfpkms$Row.names
#allfpkms <- allfpkms[-c(1)]


load(file=paste0("/wsu/home/groups/piquelab/charvey/GxE_041015/network_analysis/Full/fpkms/",platePrefixs[2],".counts.fpkm.gc.qn.adj.Rd"))
P2cv2 <- cv2
P2fpkms <- fpkms

#allfpkms <- merge(P2fpkms,allfpkms, by=0, all=FALSE)
#rownames(allfpkms) <- allfpkms$Row.names
#allfpkms <- allfpkms[-c(1)]

rm(list=c("fpkms","cv2")) #Cleanup



PBMCfpkms <- merge(P1fpkms, P2fpkms, by=0, all=FALSE)
PBMCcv2 <- smartbind(P1cv2,P2cv2)


allfpkms <- merge(P1fpkms, P2fpkms, by=0, all=FALSE)
rownames(allfpkms) <- allfpkms$Row.names
allfpkms <- allfpkms[-c(1)]


#Log 2 normalize data
allfpkms <- allfpkms/log10(2)

##Take transcripts that map to unique genes


mytP1 <- read.table(paste0("/wsu/home/groups/piquelab/gmb/GxE/differential_expression/expression/analysis/data/expr/",platePrefixs[1],".bwa.counts.fpkm.gc.topTx.txt.gz"),sep="\t",as.is=T,header=T)
rownames(mytP1) <- mytP1$t.id


## mytAll <- merge(mytP1, mytAll, by=0, all=FALSE)
## rownames(mytAll) <- mytAll$Row.names


mytP2 <- read.table(paste0("/wsu/home/groups/piquelab/gmb/GxE/differential_expression/expression/analysis/data/expr/",platePrefixs[2],".bwa.counts.fpkm.gc.topTx.txt.gz"),sep="\t",as.is=T,header=T)
 rownames(mytP2) <- mytP2$t.id

PBMCmyt <- merge(mytP1, mytP2, by=4, all=FALSE)



rownames(PBMCfpkms) <- PBMCfpkms$Row.names
PBMCfpkms <- PBMCfpkms[,!(colnames(PBMCfpkms)%in% droper)]

## mytAll <- merge(mytP2, mytAll, by=0, all=FALSE)
## rownames(mytAll) <- mytAll$Row.names



#Collapse data frame to only include 1 transcript per g.id [Gene ID}
#mytAll <- merge(mytP1, mytP2, by=4, all=FALSE)



##Read in treatment key 
treatmentKEY <- read.table("/wsu/home/groups/piquelab/gmb/GxE_full/analysis/finalFigures/treatmentKey.txt", as.is=T, sep='\t', comment.char="", header=TRUE)

#Chosing top transcript from .topTx.txt file tables
anno <- read.table("/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.bed.gz")

platePrefixs <- c("DP1","DP7")
load(file=paste0("/wsu/home/groups/piquelab/gmb/GxE_full/expression/fpkms/",platePrefixs[1],"/",platePrefixs[1],".star.counts.fpkm.gc.topTx.Rd"))
trans1 <- names(topTx)
load(file=paste0("/wsu/home/groups/piquelab/gmb/GxE_full/expression/fpkms/",platePrefixs[2],"/",platePrefixs[2],".star.counts.fpkm.gc.topTx.Rd"))
trans2 <- names(topTx)
mergetrans <- merge(trans1,trans2,by=1,all=FALSE)

indTr <- !is.na(match(rownames(LCLfpkms),mergetrans[,1]))
LCLdataExpr <- LCLfpkms[indTr,]
rownames(LCLdataExpr) <- anno$V13[match(rownames(LCLdataExpr),anno$V4)]


#HUVECdataExpr <- t(HUVECdataExpr)

####################################
indTr <- !is.na(match(rownames(HUVECfpkms),HUVECmyt$t.id))
HUVECdataExpr <- HUVECfpkms[indTr,]
rownames(HUVECdataExpr) <- HUVECmyt$ensg.x[match(rownames(HUVECdataExpr),HUVECmyt$t.id)]
HUVECdataExpr <- t(HUVECdataExpr)

#Merging the data frames by ensg, repeat for every Cell type
dataExpr <- merge(dataExpr, PBMCdataExpr, by=0, all=FALSE)
rownames(dataExpr) <- dataExpr$Row.names

dataExpr <- t(dataExpr)
#Removing Selenium From The Treatment List
#dataExpr <- dataExpr[-grep("T19C1",rownames(dataExpr)),]

#Look for outlier Samples (PDF)

sampleTree <- hclust(dist(dataExpr), method = "average");
pdf(paste0("Sample.Dendogram." ,celline, ".GC.QNorm.Adj.pdf"))
par(cex = 0.6);
par(mar = c(0,4,2,0))
Tree_treatments <- sapply(strsplit(sampleTree$labels,".",fixed=TRUE),"[",2)
labCols <- ColKey$Color[match(Tree_treatments, ColKey$ID)]
sampleTree$nodepar <- labCols
plotDendroAndColors(sampleTree, labCols, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#Remove outlier LCL dataExpr <- dataExpr[-c(9,15),] (Iron Treatment)

#Determine what soft thresholding power should be used

indTr <- !is.na(match(rownames(PBMCfpkms),mergetrans[,1]))
PBMCdataExpr <- PBMCfpkms[indTr,]
rownames(PBMCdataExpr) <- anno$V13[match(rownames(PBMCdataExpr),anno$V4)]


powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(dataExpr, powerVector=powers, verbose=5)

#Assigning a soft thresholding power: Used 12 for LCL analysis and 5 for HUVEC 12 for HUVEC without mean Adjustment, 6 for HUVEC without NaSE treatment
#Default of 6 for unsigned networks is used if there is no 

pwr <- sft$powerEstimate 
if(is.na(sft$powerEstimate))
  pwr <- 6

#Scale free topology model plot
pdf(paste0("Scale.Free.Topology.Model.",celline,".GC.QNorm.Adj.pdf"))
cex1=0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",type="n", main=paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers, cex=cex1, col="red");                                                   	
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()


#Build Network Blockwise Method (Use for large RNAseq datasets)
net <- blockwiseModules(dataExpr, maxBlockSize=40000, power = pwr, TOMType="unsigned", minModuleSize=10, numericLabels=TRUE, saveTOMs=TRUE, saveTOMFileBase= paste0(platePrefixs[1],".TOM-blockwise.GC.QNorm.Adj.topTx.TOM"), verbose=3)

#Set up network variables
MEs <- net$MEs;
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
transcripts <- colnames(dataExpr)
modVector <- c(0:max(net$colors))

#Save network data
save(list=c("dataExpr","net","pwr","sampleTree","sft","moduleColors","moduleLabels","MEs"),
	file=paste0("", celline, ".GC.QNorm.Adj.HUVEC.SignedNetwork.Rd"))

##Network Visualization

#Network Dendogram
pdf(paste0("Dendogram.",celline,".GC.QNorm.ADJ.SG.pdf"))
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
dev.off()

#Module label mappings
ModLabMap <- data.frame(unique(net$colors),unique(moduleColors))
colnames(ModLabMap) <- c("NumericLables","Colors")
tMEs <- t(MEs)
colnames(tMEs) <- rownames(dataExpr)
Individuals <- sapply(strsplit(colnames(tMEs),".",fixed=TRUE),"[",1)
Treatments <- sapply(strsplit(colnames(tMEs),".",fixed=TRUE),"[",2)
colFact <- factor(Treatments)

#Eigengene Heatmap

Rcols <- ModLabMap$Colors[match(substring(rownames(tMEs),3),ModLabMap$NumericLables)]
Ccols <- as.vector(ColKey$Color[colFact])
pdf(paste0(celline,".Eigengene.Heatmap.GC.QNorm.ADJ.SG.pdf"),width=12,height=12)
Eigengene.hmp <- heatmap.2(tMEs,trace="n",cexCol=0.9,margins=c(7,5),col=redblue, ColSideColors=Ccols)
dev.off()


#Smartbind the the covariate tables for 2 plates

allcv2 <- smartbind(P1cv2,P2cv2)

#colnames(DP8cv2) <- colnames(DP4cv2)
#HUVcv2 <- rbind(DP4cv2,DP8cv2)
#Remove the treatments that are screwing things up
#levels(LCLcv2$Treatment.ID) <- LCLcv2$Treatment.ID
#LCLcv2 <- smartbind(DP1cv2,DP7cv2)
#LCLcv2 <- LCLcv2[-c(9,15),]

BuildResponseHeatmap <- function(networkOBJ,cvTab){
  tVals <- sapply(networkOBJ$MEs,function(y){
	aux <- lm(y ~ cvTab$Plate.ID + cvTab$Treatment.ID)
	s.aux <- summary(aux)
	s.aux$coefficients[,3]
      })
  heatmapOBJ <- as.data.frame(tVals[-c(1,3,4,5),])
  return(heatmapOBJ)
}


#M1 = Plate + celltype + treatment
#M2 = plate + cell-type + treatment + cell-type X treatment
#anova(M1,M2)
tVals <- sapply(net$MEs,function(y){
	aux <- lm(y ~ allcv2$Plate.ID + allcv2$Treatment.ID + allcv2$CellType)
	s.aux <- summary(aux)
	s.aux$coefficients[,3]
      })



library(reshape2)
comb <- melt(table(allcv2$Treatment.ID,allcv2$CellType))
comb <- comb[comb$value>0,]
row.names(comb) <- paste0(comb$Var1,":",comb$Var2)

comb <- comb[!(comb$Var1 %in% c("CO1","CO2","CO3")),]

tr2ct <- unique(allcv2[,c("Treatment.ID","Control.ID")])
row.names(tr2ct)<- as.character(tr2ct$Treatment.ID)

Fnew<- sapply(net$MEs,function(y){  
        sapply(row.names(comb),function(rn){
          xt  <- allcv2$Treatment.ID == comb[rn,"Var1"] & allcv2$CellType == comb[rn,"Var2"]
          tr1 <- allcv2$Treatment.ID;
          tr0 <- tr1
          tr0[xt] <- as.character(tr2ct[comb[rn,"Var1"],"Control.ID"])
          ##& allcv2$CellType == as.comb[rn,"Var2"]
          x <- tr1!=tr0
          m0 <- lm(y ~ allcv2$Plate.ID + tr0)
          m1 <- lm(y ~ allcv2$Plate.ID + tr0 + x)
          aux <- anova(m0,m1)
          aux$F[2]
        })
      })

Zscores <- sapply(net$MEs,function(y){  
        sapply(row.names(comb),function(rn){
          xt  <- allcv2$Treatment.ID == comb[rn,"Var1"] & allcv2$CellType == comb[rn,"Var2"]
          xc  <- allcv2$Plate.ID == allcv2$Plate.ID[xt][1] & allcv2$Treatment.ID == as.character(tr2ct[comb[rn,"Var1"],"Control.ID"])
          aux <- t.test(y[xt],y[xc])
          z = -qnorm(aux$p.value/2)*sign(aux$statistic)
          as.numeric(z)
        })
      })




Zscores.m <- melt(Zscores)

colnames(Zscores.m) <- c("trAndCell","module","z")

Zscores.m <- cbind(Zscores.m, comb[Zscores.m$trAndCell,c("Var1","Var2")])


##Zscores.meanTreatment <- aggregate(Zscores.m$z,by=list(m=Zscores.m$module,tr=Zscores.m$Var1),mean)
##removeThese <- names(table(Zscores.m$Var1)[table(Zscores.m$Var1) /88 <2]
##Zscores.mNew <- Zscores.m[!Zscores.m$Var1 %in% removeThese,]
##Zscores.meanTreatment <- aggregate(Zscores.mNew$z,by=list(m=Zscores.mNew$module,tr=Zscores.mNew$Var1),mean)
##Zscores.meanTreatment <- as.matrix(acast(Zscores.meanTreatment, m ~ tr,value.var="x"))
##Zscores.meanCelltype <- aggregate(Zscores.m$z,by=list(m=Zscores.m$module,tr=Zscores.m$Var2),mean)
##Zscores.meanCelltype <- as.matrix(acast(Zscores.meanCelltype, m ~ tr,value.var="x"))


Zscores.tr <- aggregate(Zscores.m$z,by=list(m=Zscores.m$module,tr=Zscores.m$Var1),mean)

Zscores.ctype <- aggregate(Zscores.m$z,by=list(m=Zscores.m$module,ctype=Zscores.m$Var2),mean)

Zscores.tr <- as.matrix(acast(Zscores.tr, m ~ tr,value.var="x"))


whichTreatZscore.tr <- colnames(Zscores.tr)[apply(abs(Zscores.tr),1,which.max)]
names(whichTreatZscore.tr) <- sapply(strsplit(rownames(Zscores.tr),"ME"),'[',2)
ExprTbl$ZscoreCols <- col2hex(treatmentKEY[as.character(whichTreatZscore.tr[match(ExprTbl$mods,names(whichTreatZscore.tr))]),]$DeepColor)
ExprTbl$ZscoreTreatmentName <- treatmentKEY[as.character(whichTreatZscore.tr[match(ExprTbl$mods,names(whichTreatZscore.tr))]),]$Common_Name


#insignificant modules 31 and 35
numberSig <- rowSums(abs(t(Zscores))>2)[rownames(ModData)]
save(list=c("Enrichment.Heatmap","GOEnr","Zscores","numberSig","Zscores.meanCelltype","Zscores.meanTreatment"),file="../Network_Data.Rd")

##ModData[order(ModData[,2]),]

colnames(Zscores.tr) <- treatmentKEY[colnames(Zscores.tr),"Common_Name"]
Zscores.tr[(abs(Zscores.tr) < 2.57)] <- 0

pdf("Zscore.tr.Heatmap.pdf",height=14,width=14)
Fnew.HMP <- heatmap.2(Zscores.tr,trace="n",cexCol=1,,cexRow=1,margins=c(22,6),col=redblue(75), ColSideColors =treatmentKEY$DeepColor[match(colnames(Zscores.tr),treatmentKEY$Common_Name)])
dev.off()



logPvalScores <- sapply(net$MEs,function(y){  
        sapply(row.names(comb),function(rn){
          xt  <- allcv2$Treatment.ID == comb[rn,"Var1"] & allcv2$CellType == comb[rn,"Var2"]
          xc  <- allcv2$Plate.ID == allcv2$Plate.ID[xt][1] & allcv2$Treatment.ID == as.character(tr2ct[comb[rn,"Var1"],"Control.ID"])
          aux <- t.test(y[xt],y[xc])
          z = -qnorm(aux$p.value/2)*sign(aux$statistic)
          -log10(aux$p.value)*sign(aux$statistic)
        })
      })


##Zscores[abs(Zscores)<2.63] <- 0

o <- order(comb$Var1)

Zscores <- Zscores[o,]

## treatmentKEY[as.character(comb$Var1),"DeepColor"]
pdf("Zscore.Heatmap.pdf",height=14,width=14)
Fnew.HMP <- heatmap.2(t(Zscores),Colv= FALSE,trace="n",cexCol=1,,cexRow=1,margins=c(22,6),col=redblue(75), ColSideColors =treatmentKEY[as.character(comb$Var1[o]),"DeepColor"])
dev.off()


netSize <- as.numeric(table(net$colors))
hmp <- as.data.frame(tVals[-c(1),])
colnames(hmp) <- paste(gsub("ME","M",colnames(hmp))," (",netSize[(as.numeric(substring(colnames(hmp),3)) + 1)],")",sep="")
tHmp <- t(hmp)

##Build Response Heatmap

#Remove Control and Cell Plate Columns
#tHmp <- tHmp[,-c(1,2)]

colnames(tHmp) <- sapply(strsplit(colnames(tHmp),"ID",fixed=T),"[",2)

pdf(paste0("Tvalue.Heatmap.",celline,".pdf"),height=12,width=12)
Treatments <- colnames(tHmp)
colFact <- factor(Treatments)
RCols <- ModLabMap$Colors[match(substring(rownames(tHmp),3),ModLabMap$NumericLables)]
#CCols <- ColKey$Color[colFact]
lm.Heatmap <- heatmap.2(tHmp,trace="n",cexCol=0.9,margins=c(7,5),col=redblue(75))
dev.off()


#####
pdf("Tvalue.Heatmap.Final.pdf",height=12,width=12)
Treatments <- colnames(FinalHMP)
colFact <- factor(Treatments)
RColgs <- ModLabMap$Colors[match(substring(rownames(FinalHMP),3),ModLabMap$NumericLables)] 
lm.Heatmap <- heatmap.2(FinalHMP,trace="n",cexCol=0.9,margins=c(7,5),col=redblue(75))
dev.off()
###


#Build adjacency matrix and set up ExprTbl data frame
AdjM <- unsignedAdjacency(dataExpr, power=pwr)

#load(net$TOMFiles)
#asM.TOM <- as.matrix(TOM)
Alldegrees1 <- intramodularConnectivity(AdjM, moduleColors)
ExprTbl <- Alldegrees1
ExprTbl$g.id <- mytP1$g.id[match(colnames(dataExpr),mytP1$ensg)]
ExprTbl$mods <- net$colors
ExprTbl$modCol <- moduleColors
rownames(ExprTbl) <- colnames(dataExpr)
ExprTbl$t.id <- anno$V4[match(rownames(ExprTbl),anno$V13)]DEG
#ExprTbl$ensg <- mytAll$ensg.x[match(rownames(ExprTbl),mytAll$t.id)]

##Eigengene Dissimilarity Tree 
#dissimME=(1-t(cor(MEs, method="p")))/2
#hclustdatME=hclust(as.dist(dissimME), method="average" )
#datMEcolors <- sapply(hclustdatME$labels, function(xx){
#  modnum <- substring(xx,3)
#    modtreat <- ExprTbl$Treatment[match(modnum,ExprTbl$mods)]
#  cat(modtreat)
#  ColIDX <-match(substring((strsplit(modtreat,".",fixed=T)[[1]][2]),3),ColKey$ID)#  ColKey$Color[ColIDX]
#})


##Gene Significance
## REMOVED FROM EXPRTBL
#Determine Gene significance by fitting to linear model
#GenetVals <- ParallelSapply(dataExpr,function(y){
#	aux <- lm(y ~ HUVcv2$Plate.ID + HUVcv2$CellLine + HUVcv2$Treatment.ID)
#	s.aux <- summary(aux)
#	s.aux$coefficients[,3]
#})
#tValsPergene <- data.frame(GenetVals)
#Now add Treatment and Gene Significance Column to ExprtTbl Data frame

ModResponseTval <- rownames(hmp[-c(1:9),])[apply(abs(hmp[-c(1:9),]), 2, which.max)]
names(ModResponseTval) <- sapply(strsplit(sapply(strsplit(colnames(hmp),"M"),'[',2)," "),'[',1)
ExprTbl$TvalTreatment <- sapply(strsplit(ModResponseTval[match(ExprTbl$mods,names(ModResponseTval))],"ID"),'[',2)

#TreatIDX <- match(ExprTbl$Treatment, rownames(tValsPergene))
#ExprTbl$GeneSig <- abs(t(tValsPergene)[cbind(seq_along(TreatIDX),TreatIDX)])
## ModTreatmentResponse <- rownames(log2Enrichment)[apply(abs(log2Enrichment),2,which.max)]
## names(ModTreatmentResponse) <- colnames(log2Enrichment)
## ExprTbl$Treatment <- as.character(ModTreatmentResponse[match(ExprTbl$mods,substring(names(ModTreatmentResponse),2))])

## #ModData data frame [A Data Frame with information about network modules]
## ModData <- as.data.frame(t(log2Enrichment))

#####
##ModData$Specificity <- apply(ModData,1,function(y) {
##	ddd <- max(y)/length(y)
##})
#####

ModData <- data.frame(t(hmp))
ModData$Specificity <- apply(ModData,1,function(y) {
	ddd <- max(abs(y))/sum(abs(y))
})

colnamer <- sapply(strsplit(colnames(ModData),".",fixed=T),function(ll){
  if (length(ll) > 1){
    substring(ll[3],3)
  }
  else {
    ll[1]
  }
})
colnamer[1] <- "PlateEffect"
colnames(ModData) <- colnamer


##Gene significance per DEG data
DEGdf <- as.data.frame(ExprTbl$TvalTreatment)
colnames(DEGdf) <- "Treatment"
rownames(DEGdf) <- rownames(ExprTbl)

SuperDEdf <- 

allDEfiles <- lapply(plateKey, function(lister){
  filelist <- list.files(paste0("/wsu/home/groups/piquelab/gmb/GxE_full/expression/DESeq2/out_data_",lister,"/stats"),full.names=TRUE)
  names(filelist) <- gsub("\\.txt","",gsub(".*_","",filelist))
  filelist
})
  sapply(filelist, function(xx){
    cat(paste0(xx,"\n"))
    final <- read.table(xx,header=T,as.is=T,sep=" ",row.names="t.id")
    as.data.frame(final)
  })
})

filelist <- as.character()
for(ff in platePrefixs){
  filelist <- c(filelist,list.files(paste0("/wsu/home/groups/piquelab/gmb/GxE_full/expression/DESeq2/out_data_",ff,"/stats"),full.names=T))
}

for(qq in filelist){
  tabler <- read.table(qq,header=T)
  tname <- strsplit(strsplit(qq,"_")[[1]][8],".",fixed=T)[[1]][1]
  newcol <- -log10(tabler$padj[match(rownames(DEGdf),tabler$t.id)])
  DEGdf[,paste(tname)] <- newcol
}

names(filelist) <- gsub("\\.txt","",gsub(".*_","",filelist))
filelist <- filelist[names(filelist) != "CO2"]
filelist <- filelist[names(filelist) != "CO1"]

DEdfs <- lapply(cellines,function(cellineItter){
   cat(cellineItter,"\n")
   lapply(cellineItter, function(filelist)
   tabler <- read.table(filer,header=T,as.is=T,sep=" ",row.names="t.id")
   tabler <- tabler[rownames(ExprTbl)%in% tabler$ensg,]
   tabler$mods <- ExprTbl[,"mods"]
   tabler
 })
names(modVector) <- paste0("M",modVector)

BinResults <- lapply(modVector,function(mod){
  lapply(DEdfs,function(tabler){
    pvals <- tabler[tabler$mods==mod,"pval"]
    npos <- sum(pvals<=0.01,na.rm=T)+1
    nneg <- sum(pvals>0.01,na.rm=T)+1
    bt <- binom.test(npos,nneg+npos,0.01)
    bt
  })
})


FishResults <- lapply(modVector,function(mod){
  lapply(DEdfs,function(tabler){
    tt <- table(tabler$mods==mod,tabler[,"padj"]<0.1)+1
    bt <- fisher.test(tt)
    bt
  })
})


log2Enrichment <- sapply(FishResults,function(xx){
  sapply(xx,function(bt){
    log2(as.numeric(bt$estimate))
  })
})


##Add Fisher-identified Treatment to ExprTbl

whichTreatZscore <- sapply(strsplit(rownames(Zscores)[apply(abs(Zscores),2,which.max)],":"),"[",1)
names(whichTreatZscore) <- sapply(strsplit(colnames(Zscores),"ME"),'[',2)
ExprTbl$ZscoreCols <- col2hex(treatmentKEY[as.character(whichTreatZscore[match(ExprTbl$mods,names(whichTreatZscore))]),]$DeepColor)


whichTreatFish <- rownames(log2Enrichment)[apply(log2Enrichment,2,which.max)]
names(whichTreatFish) <- sapply(strsplit(colnames(log2Enrichment),"M"),'[',2)
ExprTbl$FishTreatment <- whichTreatFish[match(ExprTbl$mods,names(whichTreatFish))]

##################################################################################
## htestColnames <- names(BinResults[[1]][[1]])                                 ##
## ListParser <- lapply(BinResults,function(xx){                                ##
##   outer <- as.data.frame(sapply(xx,cbind))                                   ##
##   rownames(outer) <- htestColnames                                           ##
##   outer                                                                      ##
## })                                                                           ##
##                                                                              ##
## Estimates <- as.data.frame(sapply(ListParser, function(yy){                  ##
##   estimate <- unlist(yy[5,])                                                 ##
## }))                                                                          ##
##                                                                              ##
##                                                                              ##
## colnames(Estimates) <- paste0("M",modVector)                                 ##
## rownames(Estimates) <- as.character(sapply(rownames(Estimates),function(zz){ ##
##   strsplit(zz,".",fixed=T)[[1]][1]                                           ##
## }))                                                                          ##
##################################################################################

pdf("Tvalue.Heatmap.pdf",height=12,width=12)
Treatments <- substring(sapply(strsplit(colnames(tHmp),".",fixed=TRUE),"[",2),3)
colFact <- factor(Treatments)
RCols <- ModLabMap$Colors[match(substring(rownames(tHmp),3),ModLabMap$NumericLables)]
CCols <- ColKey$Color[colFact]
response.Heatmap <- heatmap.2(tHmp,trace="n",cexCol=0.9,margins=c(7,5),col=redblue(75),ColSideColors=CCols)
dev.off()



pdf("Fisher.Test.Heatmap.pdf",height=12,width=12)
Fisher.Heatmap <- heatmap.2(t(log2Enrichment),trace="n",cexCol=0.9,cexRow=0.6,margins=c(7,5),col=redblue)
dev.off()

##Reordering the rows
## pdf("Fisher.Test.Heatmap.pdf",height=12,width=12)
## ReorderColumns <- as.numeric(substring(colnames(MEs),3)) + 1
## Rowidx <- response.Heatmap$rowInd
## Colidx <- response.Heatmap$colInd
## HmpMat <- t(as.matrix(log2Enrichment[colnames(tHmp),ReorderColumns])[Colidx,rev(Rowidx)])
## Fisher.Heatmap <- heatmap.2(HmpMat,trace="n",Rowv=F,Colv=F,dendrogram="none",cexCol=0.9,cexRow=0.6,margins=c(7,5),col=redblue)
## dev.off()


#### Block builds a heatmap comparison
## pdf("~/share/HUVEC/NO_NaSE/Compare.heatmap.pdf")
## layout(rbind(c(4,3,8,7),c(2,1,6,5)),
##               widths = c(1,2,1,2), heights = c(1,2), respect = FALSE)
## response.Heatmap <- heatmap.multi(tHmp,trace="n",dendrogram="none",cexRow=0.6,cexCol=0.9,margins=c(7,5),col=redblue(75))
## Fisher.Heatmap <- heatmap.multi(HmpMat,trace="n",Rowv=F,Colv=F,dendrogram="none",cexCol=0.9,cexRow=0.6,margins=c(7,5),col=redblue)
## dev.off()



#Plotting out the Connectivity vs Differential Expression 
## numMods <- length(unique(moduleLabels))

## pdf(paste0(celline,"TOMconnectivity.Vs.Significance.pdf"))

## for(mm in modVector){
##   modIDX  <- net$colors == mm
##   responseTreat <- ExprTbl$TvalTreatment[modIDX][1]
##   modcol <- moduleColors[modIDX][1]
##   plot(as.numeric(DEGdf[modIDX,responseTreat])+1e-144,as.numeric(ExprTbl$kWithin[modIDX]),col=modcol,main=paste0(mm,"Hubby Vs Intramodular Connectivity"),xlab="-log(pval Adjusted)",ylab="Intramodular Connectivity")
##   text(ExprTbl$logPadj[modIDX],ExprTbl$kWithin[modIDX],ExprTbl$g.id[modIDX],col="black",cex=0.4)
## }
## dev.off()


###########################################################################################
## #This part can be run interactively, choising the modules and treatments of interest
## treatment <- "T6C1"     #Loads DEX data for the treatment of choice            ##
## OfInterestMods <- (23)  #Chose which module you would like to look at          ##
## threshold <- 0.1       #Chose the threshold to use                             ##
###########################################################################################

#Calculate datKME to give module significance measure for each gene
dataKME <- signedKME(dataExpr,net$MEs,outputColumnName="KK.")
kkIDX <- match(ExprTbl$mods, substring(colnames(dataKME),4))
ExprTbl$corME <- dataKME[cbind(seq_along(kkIDX),kkIDX)]

## #Build All Cytoscape Modules Using Fisher Test Treatments
maxT <- apply(abs(hmp),2,max) 
names(maxT) <- substring(sapply(strsplit(names(maxT), " "),"[",1),2)
OfInterestMods <- as.numeric(substring(sapply(strsplit(colnames(hmp)[apply(abs(hmp),2,max) > 8], " "), "[", 1),2))

ExprTbl$TreatmentColor <- ColKey$Color[match(ExprTbl$TvalTreatment,ColKey$ID)]
ExprTbl$maxT <- as.character(maxT)[match(ExprTbl$mods,names(maxT))]

system("mkdir CytoscapeFiles")
setwd("./CytoscapeFiles")
threshold <- 0.05
Cytoscape <- sapply(as.numeric(names(whichTreatFish)), function(OfInterestMods){
  indexer <- OfInterestMods + 1
  treatment <- whichTreat[indexer]
  result <- c(OfInterestMods, treatment, threshold)
    cyto <- Vizmodule(OfInterestMods, treatment, threshold, celltype)
})


lapply
setwd("../")
rm(Cytoscape)


#Color annotations for network figure
sapply(abs(hmp), function(fff){
  if max(fff) > 8{
    

#Build All Cytoscape Modules Using T-value Treatments
OfInterestMods <- modVector
ExprTbl$test <- as.numeric(hmp[18,][match(ExprTbl$mods,substring(sapply(strsplit(names(hmp[18,])," "),'[',1),2))])
threshold <- 0.05
treatment <- "T6C1"
system("mkdir CytoscapeFiles")
setwd("./CytoscapeFiles")
threshold <- 0.05
Cytoscape <- sapply(as.numeric(sort(unique(ExprTbl$mods))), function(OfInterestMods){
  treatment <- ExprTbl$TvalTreatment[ExprTbl$mods == OfInterestMods][1]
     cyto <- Vizmodule(OfInterestMods, treatment, threshold)
})
setwd("../")
rm(Cytoscape)


#Change the module assignment by using a dynamic tree cutting method

#TreeCutter <- cutreeDynamicTree(net$dendrograms[[1]],minModuleSize=15)
#TreeCutterColors <- labels2colors(TreeCutter)
#TreeCutterEigengenes <- moduleEigengenes(dataExpr,TreeCutter, softPower=5)

#pdf("~/share/HUVEC/TreeCut/Average.Expression.Heatmap.pdf",width=12,height=12) #View Heatmap of Average Module Expression
#heatmap.2(t(TreeCutterEigengenes$averageExpr),trace="n",cexCol=0.9,cexRow=0.6,margins=c(7,5),col=redblue)
#dev.off()

#Set up data frame for significant module heatmap
## moduleSort <- names(sort(apply(SortModData, 1, max),decreasing=TRUE))[1:20]
## moduleSort <- substring(moduleSort, 3)
## GeneHeatmap <- ExprTbl[ExprTbl$mods %in% moduleSort,]
## OrderedGeneHeatmap <- GeneHeatmap[order(GeneHeatmap$mods),]

#treecutter3 <- cutreeDynamic(net$dendrograms[[1]],cutHeight=0.995, method='hybrid',minClusterSize=10,distM=matTOM,deepSplit=3)
#mod3Cols <- lables2colors(treecutter3)
#MEs3 <- moduleEigengenes(dataExpr,mod3Cols,softPower=6)
#Enrichment Analysis of Modules in the network

mart <- useDataset("hsapiens_gene_ensembl", useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL"))
GeneData <- getBM(filters="ensembl_gene_id",attributes= c("ensembl_gene_id","entrezgene"),values=rownames(ExprTbl),mart=mart)

EntrezgeneIDs <- GeneData$entrezgene[match(rownames(ExprTbl),(unique(GeneData$ensembl_gene_id)))]
GOEnr <- GOenrichmentAnalysis(labels=net$colors,entrezCodes=EntrezgeneIDs,organism="human",ontologies=c("BP"),backgroundType="allGiven",verbose=1)
EnrichData <- as.data.frame(rbind(GOEnr$bestPTerms[[1]]$enrichment))
EnrichData$negLOGp <- -log10(EnrichData$enrichmentP)



cc <- heatmap.3(lm.Heatmap)

namesOfInterest <- unique(EnrichData$termID[EnrichData$negLOGp > 3.5])
test <- as.data.frame(GOEnr$enrichmentP)
Enrichment.Heatmap <- test[,match(namesOfInterest,colnames(test))]
Enrichment.Heatmap <- -log10(Enrichment.Heatmap)
Enrichment.Names <- unique(EnrichData$termName[EnrichData$negLOGp > 3.5])
colnames(Enrichment.Heatmap) <- Enrichment.Names
Enrichment.Heatmap <- Enrichment.Heatmap[substring(sapply(strsplit(rownames(tHmp)," "),"[",1),2),]
colmaker <- colorRampPalette(c("white","red"))
rownames(Enrichment.Heatmap) <- rownames(Final.HMP)
Enrichment.Heatmap <- Enrichment.Heatmap[cc$rowInd,]
Enrichment.Heatmap <- Enrichment.Heatmap[c(88:1),]

pdf('GOheatmap2.pdf', height=14, width=14)
GOhmp <- heatmap.2(dendrogram = "column",as.matrix(Enrichment.Heatmap),Rowv=FALSE,trace="n",cexCol=0.5,cexRow=1,margins=c(22,6),col=colmaker(15),RowSideColors=rowCols)
dev.off()



########


nTrRep <- table(allcv2$Treatment.ID)
trNames <- names(nTrRep[nTrRep>0 & nTrRep <30])
names(trNames) <- trNames

##Association with module eigengene with treatment that takes into account interaction with cell type
pTreatAll<- sapply(net$MEs,function(y){
        sapply(trNames, function(rn){
          tr <- allcv2$Treatment.ID;
          x  <- tr == rn;
          tr[x] <- "CO1"
          m0 <- lm(y ~ allcv2$Plate.ID + tr)
          m1 <- lm(y ~ allcv2$Plate.ID + tr + allcv2$CellType * x)
          aux <- anova(m0,m1)
##          aux$F[2]
          -log10(aux$"Pr(>F)"[2])
        })
      })

F.TreatAll<- sapply(net$MEs,function(y){
        sapply(trNames, function(rn){
          tr <- allcv2$Treatment.ID;
          x  <- tr == rn;
          tr[x] <- "CO1"
          m0 <- lm(y ~ allcv2$Plate.ID + tr)
          m1 <- lm(y ~ allcv2$Plate.ID + tr + allcv2$CellType * x)
          aux <- anova(m0,m1)
          aux$F[2]
          ##-log10(aux$"Pr(>F)"[2])
        })
      })




rownames(pTreatAll) <- treatmentKEY$Common_Name[match(rownames(pTreatAll),rownames(treatmentKEY))]
pTreat.colcolors <- treatmentKEY$DeepColor[match(rownames(pTreatAll),treatmentKEY$Common_Name)]
pdf("Treatment.Heatmap.pdf",height=14,width=14)
pTreat.HMP <- heatmap.2(t(pTreatAll),trace="n",cexCol=2,cexRow=1,margins=c(22,6),col=colmaker(75),ColSideColors=pTreat.colcolors)
dev.off()




##cell-type x Treatment
nTrRep <- table(allcv2$Treatment.ID)
trNames <- names(nTrRep[nTrRep>3 & nTrRep <30])
names(trNames) <- trNames
pInter<- sapply(net$MEs,function(y){
        m0 <- lm(y ~ allcv2$Plate.ID + allcv2$Treatment.ID)
        sapply(trNames, function(rn){
          x  <- allcv2$Treatment.ID == rn;
          m1 <- lm(y ~ allcv2$Plate.ID + allcv2$Treatment.ID + allcv2$CellType * x)
          aux <- anova(m0,m1)
##          aux$F[2]
          -log10(aux$"Pr(>F)"[2])
        })
      })

rownames(pInter) <- treatmentKEY$Common_Name[match(rownames(pInter),rownames(treatmentKEY))]
pInter.colcolors <- treatmentKEY$DeepColor[match(rownames(pInter),treatmentKEY$Common_Name)]
pdf("Celltype_Treatment_interaction.Heatmap.pdf",height=14,width=14)
pInter.HMP <- heatmap.2(t(pInter),trace="n",cexCol=2,cexRow=1,margins=c(22,6),col=colmaker(75),ColSideColors=pInter.colcolors)
dev.off()


####


##Deep split analysis

## netDeep1 <- blockwiseModules(dataExpr, maxBlockSize=40000, power = 6, TOMType="unsigned", minModuleSize=10, deepSplit=1, numericLabels=TRUE, saveTOMs=TRUE, saveTOMFileBase= paste0("",platePrefix,".TOM-blockwise.power12.GC.QNorm.Adj.SingleGene.HUVEC"), verbose=3)
## netDeep3 <- blockwiseModules(dataExpr, maxBlockSize=40000, power = 6, TOMType="unsigned", minModuleSize=10, deepSplit=3, numericLabels=TRUE, saveTOMs=TRUE, saveTOMFileBase= paste0("",platePrefix,".TOM-blockwise.power12.GC.QNorm.Adj.SingleGene.HUVEC"), verbose=3)
## netDeep4 <- blockwiseModules(dataExpr, maxBlockSize=40000, power = 6, TOMType="unsigned", minModuleSize=10, deepSplit=4, numericLabels=TRUE, saveTOMs=TRUE, saveTOMFileBase= paste0("",platePrefix,".TOM-blockwise.power12.GC.QNorm.Adj.SingleGene.HUVEC"), verbose=3)

## GetSpecificity <- function(netObj,cvTabel){
##     tValstmp <- sapply(netObj$MEs,function(y){
## 	aux <- lm(y ~ HUVcv2$Plate.ID + HUVcv2$CellLine + HUVcv2$Treatment.ID)
## 	s.aux <- summary(aux)
## 	s.aux$coefficients[,3]
##     })
##     fff <- as.data.frame(tValstmp[-c(1,3,4,5),])
##     ggg <- data.frame(t(fff))
##     Specificity <- apply(ggg,1,function(y) {
##         hhh <- max(abs(y))/sum(abs(y))
##     })
##     return(Specificity)
## }

## ##Build Specificity Histogram

## Spec1 <- GetSpecificity(netDeep1,HUVcv2)
## Spec2 <- GetSpecificity(net,HUVcv2)
## Spec3 <- GetSpecificity(netDeep3,HUVcv2)
## Spec4 <- GetSpecificity(netDeep4,HUVcv2)
## hisXlim <- c(0,max(c(Spec1,Spec2,Spec3,Spec4)))
## pdf(paste0(platePrefix,".Specificity.Hisogram.pdf"))
## hist(Spec1, col="Black", xlim=hisYlim, ylim=c(0,100), breaks=50, main="Specificity Per Module Cutting Method")
## hist(Spec2, col="Red", add=T)
## hist(Spec3, col="Blue", add=T)
## hist(Spec4, col="Green", add=T)
## dev.off()


##Response Modules VS DEGs Fisher test
tId_motif <- read.table("/wsu/home/groups/piquelab/gmb/genomicAnnotation/proxGenes/huvec.tss.tfbs.txt",header=T)

mot_Tab <- table(tId_motif$pwm)

mot_Tab <- mot_Tab[mot_Tab>200]

##Mot_vector <- unique(tId_motif$pwm)
Mot_vector <- names(mot_Tab)

TfEnrich <- lapply(modVector, function(ModItter) {
     aux <-  mclapply(Mot_vector, mc.cores = 11,  function(MotItter){
                Tid_In_Module <- rownames(ExprTbl)[ExprTbl$mods == ModItter]
                        Tid_NotIn_Module <- rownames(ExprTbl)[!ExprTbl$mods == ModItter]
                        AllTid_In_Mot <- tId_motif$t.id[tId_motif$pwm == MotItter]
                        ##AllTid_NotIn_Mot <- tId_motif$t.id[!tId_motif$pwm == MotItter]
                        aa <- sum(Tid_In_Module %in% AllTid_In_Mot) + 1     #           Motif within 5KB
                        bb <- sum(Tid_NotIn_Module %in% AllTid_In_Mot) + 1 #               Y   N
                        cc <- sum(!Tid_In_Module %in% AllTid_In_Mot) + 1    #   In      Y   aa  cc
                        dd <- sum(!Tid_NotIn_Module %in% AllTid_In_Mot) + 1 #   Module  N   bb  dd
                        fisher_mat = rbind(c(aa,cc),c(bb,dd))
                        fResult <- fisher.test(fisher_mat)
                        fResult$module <- ModItter
                        fResult$pwm <- MotItter
                        fResult$fisher_mat <- fisher_mat
                        fResult
              })
      names(aux)<- Mot_vector
      aux
    })
names(TfEnrich) <- paste0("M",modVector)


##tId_motif$mods <- ExprTbl[tId_motif$t.id,"mods"]
##myTble <- table(tId_motif$mods,tId_motif$pwm)

##modNum <- table(tId_motif$mods)
##pwmNum <- table(tId_motif$pwm)

##modNum <- rowSums(myTable)
##pwmNum <- colSums(myTable)


TfEnrich.log2Odds <- sapply(TfEnrich,function(x){
  sapply(x,function(y){
    log2(y$estimate)
  })
})

TfEnrich.pVal <- sapply(TfEnrich,function(x){
  sapply(x,function(y){
    y$p.value
  })
})

#Treatment_Order <- order(as.numeric(substring(names(datMEcolors),3)),modVector)
#heatmapClmcols <- datMEcolors[Treatment_Order]
pdf("TF.enrichment.heatmap.pdf",width=12,height=12)
heatmap.2(TfEnrich.log2Odds,trace="n",cexCol=0.9,margins=c(7,5),col=redblue(75),cexRow=0.0001)
dev.off()


colnames(hmp) <- gsub("ME","M",colnames(hmp))
hmp <- hmp[,colnames(TfEnrich.log2Odds)] 



save(list=c("dataExpr","ExprTbl","BinResults","allcv2","FishResults","moduleColors","celline","DEGdf","EnrichData","log2Enrichment","net","pwr","ModData"),file=paste0(celline,".Network.Rd"))

write.table(ExprTbl, file=paste0(celline,".Network.Table.txt"))
write.table(ModData, file=paste0(celline,".Module.Table.txt"))
write.table(EnrichData, file=paste0(celline,".GO.Enrichment.txt"))
write.table(DEGdf, file=paste0(celline,".DE.Enrichment.txt"))

####
## pdf("~/share/HUVEC/NO_NaSE/Stacked.heatmap.pdf")
## h2 <- heatmap.2(TfEnrich.log2Odds,trace="n",cexCol=0.2,margins=c(7,5),col=redblue(75),cexRow=0.001)
## h3 <- heatmap.2(as.matrix(hmp),Colv=h2$colDendrogram,trace="n",cexCol=0.2,margins=c(7,5),col=redblue(75),dendrogram='none')
## sizeVec <- as.numeric(table(net$colors))
## ##layout(matrix(c(1,2,3),ncol=1),heights=c(1/2,1/4,1/4))
## ##image(TfEnrich.log2Odds[,h2$colInd],col=redblue(75))
## ##plot(sizeVec[h2$colInd],type='h',axes=F)
## dev.off()
## ###


## pdf("~/share/HUVEC/NO_NaSE/Stacked.heatmap.pdf")
## layout(rbind(c(4,3,8,7),c(2,1,6,5)),
##               widths = c(1,2,1,2), heights = c(1,2), respect = FALSE)
## test <- heatmap.multi(t(as.matrix(Estimates)),dendrogram="both",trace="n",cexCol=0.9,cexRow=0.5,margins=c(7,5),col=redblue)
## heatmap.multi(t(TfEnrich.log2Odds),Rowv=test$rowInd,trace="n",cexRow=0.5,margins=c(7,5),col=redblue(75),cexCol=0.0001)
## dev.off()
