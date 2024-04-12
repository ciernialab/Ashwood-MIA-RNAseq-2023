
#DEG analysis of RNAseq
#AVC 9/23/19

############################################################################################################################################

## code for taking counts.txt file outputs from FeatureCounts and importing into R for EdgeR and/or limmavoom analysis
## loops through each txt file and takes the GeneID and accepted_hits.bam counts columns (counts per gene)
## writes a new text file samplename.out.txt that contains only these two columns
## these files are then read into EdgeR 
## a targets dataframe is made that contains a list of the .out.txt files and the group they belong to
## the targets dataframe is then fed into EdgeR using readDGE function of EdgeR


library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(tidyr)
library(xlsx)
#loop for combining GeneID and counts for each txt file
path ="/Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts"
setwd(path)

############################################################################################################################################
#counts without UMI collapsed
############################################################################################################################################


#get file names for all samples for only UMI collapsed:
file.names <- dir(path, pattern =".collapsedUMI.counts.txt$")
file.names

#get genes
genes <- read.table(file.names[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
out.file <- as.data.frame(genes$Geneid)
colnames(out.file) <- c("Geneid")


for(i in 1:length(file.names)){
  file <- read.table(file.names[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
  counts <- data.frame(file$Geneid,file[,7])
  #extract file name without count.txt
  m <- as.data.frame(file.names[i])
  names(m) <-  c("split")
  m$split <- as.character(m$split)
  m <-tidyr::separate(m,split,into= c("name","stuff","stuff2"),sep="\\.")
  m$name <- gsub("_","",m$name)
  
  print(m$name)
  #replace column headings with names
  names(counts) <- c("Geneid",m$name)
  
  #write to output file
  out.file <- merge(out.file, counts, by= "Geneid")
}

#write out raw counts matrix:
setwd("/Users/aciernia/Sync/collaborations/Ashwood/MIA_placenta_brainRNAseq/Analysis/DEGs_fullmodel_uterinehornplus")

write.csv(out.file,"AllDataUMIcollapse.countmatrix.csv")

##The following code is for analysis of RNASeq in Limma Voom package
##########################################################################################

#Read in experient information sheet
info <- read.csv("MIARNAseq_sampleinfo.csv")


#make counts into matrix
countmatrix <- out.file[,3:ncol(out.file)]
rownames(countmatrix) <- out.file$Geneid

#order info to match count matrix
colorder <- info$SampleID
colorder <- as.character(colorder)
colorder
colnames(countmatrix)
countmatrix <- countmatrix[,colorder]
colnames(countmatrix)

countmatrix <- as.matrix(countmatrix)

##########################################################################################
#make DEG list object in EdgeR
##########################################################################################
# Make DEGList
library(edgeR)
library(limma)

#define group
info$group <- paste(info$treatment,info$Tissue, sep="_")

RG <- DGEList(counts = countmatrix, group=info$group)
# number of unique transcripts
dim(RG)
#53801    64

table(info$group)
#Poly IC_Brain Poly IC_Placenta     Saline_Brain  Saline_Placenta 
#16               16               16               16


###########################################################
# Filter for one count per million in at least 16 libraries
keep <- rowSums(cpm(RG)>1)>=16
RGcpm <- RG[keep,]
dim(RGcpm)
# 16546    64

table(rowSums(RGcpm$counts==0)==16)
geneid <- rownames(RGcpm) #ensemble IDs for biomart

#save
save(RG,RGcpm,file="DEGlist.RData")
##################################################################################
#read in gene info
##################################################################################

load("mm10_Ensemble_geneInfo.RData")

#match order: df[match(target, df$name),]
genes <- genes[match(geneid,genes$ensembl_gene_id),]


##########################################################################
#add gene info to DEGlist object
##########################################################################

#add into DEGlist
RGcpm$genes <- genes

#save
save(RG,RGcpm,file="DEGlist.RData")

###########################################################
#filtering plot
###########################################################
library(RColorBrewer)
nsamples <- ncol(RGcpm)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)


#plot:
pdf('FilteringCPM_plots.pdf', h=6,w=10)
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(RG, log=TRUE, prior.count=2)


plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")

#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(RGcpm, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()



##################################################################################
#MDS plots for all samples:
##################################################################################
info$Tissue <- as.factor(info$Tissue)
#by Tissue:
pdf(file = "MDSplot_Tissue.pdf", wi = 12, he = 10, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for treatment
levels(info$Tissue)
col.cell <- c("purple","green")[info$Tissue]
data.frame(info$Tissue,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(v,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("purple","green"),
       legend=levels(info$Tissue),
       cex = 0.8)
# Add a title
title("Tissue MDS Plot")
dev.off()


pdf(file = "MDSplot_Treatment.pdf", wi = 12, he = 10, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
info$treatment <- as.factor(info$treatment)
levels(info$treatment)
col.cell <- c("orange","skyblue")[info$treatment]
data.frame(info$treatment,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(v,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("orange","skyblue"),
       legend=levels(info$treatment),
       cex = 0.8)
# Add a title
title("Treatment MDS Plot")
dev.off()


##################################################################################
#MDS plots for brain only:
##################################################################################
#Subset Elist output from voom for only brain samples
#expr[ expr$genes$SYMBOL %in% c("Gene2", "Gene4"), ]
brainsamplesinfo <- info %>% filter(Tissue == "Brain")
brainsamples <- info %>% filter(Tissue == "Brain") %>% dplyr::select(SampleID)
brainsamples$SampleID <- as.character(brainsamples$SampleID)
exprbrain <- v[,colnames(v$E) %in% brainsamples$SampleID]

pdf(file = "MDSplot_Brain.Treatment.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(brainsamplesinfo$treatment)
col.cell <- c("orange","skyblue")[brainsamplesinfo$treatment]
data.frame(brainsamplesinfo$treatment,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprbrain,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("orange","skyblue"),
       legend=levels(brainsamplesinfo$treatment),
       cex = 0.8)
# Add a title
title("Treatment MDS Plot")
dev.off()


pdf(file = "MDSplot_Brain.sex.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(brainsamplesinfo$sex)
col.cell <- c("pink","blue")[brainsamplesinfo$sex]
data.frame(brainsamplesinfo$sex,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprbrain,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("pink","blue"),
       legend=levels(brainsamplesinfo$sex),
       cex = 0.8)
# Add a title
title("Sex MDS Plot")
dev.off()

pdf(file = "MDSplot_Brain.uterine.horn.position .pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(brainsamplesinfo$uterine.horn.position )
col.cell <- c("pink","blue","darkgreen","purple","orange","skyblue")[brainsamplesinfo$uterine.horn.position ]
data.frame(brainsamplesinfo$uterine.horn.position ,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprbrain,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill= c("pink","blue","darkgreen","purple","orange","skyblue"),
       legend=levels(brainsamplesinfo$uterine.horn.position),
       cex = 0.8)
# Add a title
title("Uterine Horn Position MDS Plot")
dev.off()

pdf(file = "MDSplot_Brain.dam.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(brainsamplesinfo$Dam)
col.cell <- getPalette(nlevels(brainsamplesinfo$Dam))[brainsamplesinfo$Dam]
data.frame(brainsamplesinfo$uterine.horn.position ,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprbrain,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill= getPalette(nlevels(brainsamplesinfo$Dam)),
       legend=levels(brainsamplesinfo$Dam),
       cex = 0.8)
# Add a title
title("Dam ID MDS Plot")
dev.off()

##################################################################################
#MDS plots for Placenta only:
##################################################################################
#Subset Elist output from voom for only Placenta samples
#expr[ expr$genes$SYMBOL %in% c("Gene2", "Gene4"), ]
Placentasamplesinfo <- info %>% filter(Tissue == "Placenta")
Placentasamples <- info %>% filter(Tissue == "Placenta") %>% dplyr::select(SampleID)
Placentasamples$SampleID <- as.character(Placentasamples$SampleID)
exprPlacenta <- v[,colnames(v$E) %in% Placentasamples$SampleID]

pdf(file = "MDSplot_Placenta.Treatment.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(Placentasamplesinfo$treatment)
col.cell <- c("orange","skyblue")[Placentasamplesinfo$treatment]
data.frame(Placentasamplesinfo$treatment,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprPlacenta,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("orange","skyblue"),
       legend=levels(Placentasamplesinfo$treatment),
       cex = 0.8)
# Add a title
title("Treatment MDS Plot")
dev.off()


pdf(file = "MDSplot_Placenta.sex.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(Placentasamplesinfo$sex)
col.cell <- c("pink","blue")[Placentasamplesinfo$sex]
data.frame(Placentasamplesinfo$sex,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprPlacenta,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill=c("pink","blue"),
       legend=levels(Placentasamplesinfo$sex),
       cex = 0.8)
# Add a title
title("Sex MDS Plot")
dev.off()

pdf(file = "MDSplot_Placenta.uterine.horn.position .pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(Placentasamplesinfo$uterine.horn.position )
col.cell <- c("pink","blue","darkgreen","purple","orange","skyblue")[Placentasamplesinfo$uterine.horn.position ]
data.frame(Placentasamplesinfo$uterine.horn.position ,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprPlacenta,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill= c("pink","blue","darkgreen","purple","orange","skyblue"),
       legend=levels(Placentasamplesinfo$uterine.horn.position),
       cex = 0.8)
# Add a title
title("Uterine Horn Position MDS Plot")
dev.off()

pdf(file = "MDSplot_Placenta.dam.pdf", wi = 8, he = 6, useDingbats=F)
#par(mfrow=c(1,2))
#plot MDS for condition
levels(Placentasamplesinfo$Dam)
col.cell <- getPalette(nlevels(Placentasamplesinfo$Dam))[Placentasamplesinfo$Dam]
data.frame(Placentasamplesinfo$uterine.horn.position ,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(exprPlacenta,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(1,1.5,
       fill= getPalette(nlevels(Placentasamplesinfo$Dam)),
       legend=levels(Placentasamplesinfo$Dam),
       cex = 0.8)
# Add a title
title("Dam ID MDS Plot")
dev.off()

###########################################################
###########################################################
#reset library sizes
RGcpm$samples$lib.size <- colSums(RGcpm$counts)


#plot library sizes
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(RGcpm$samples$lib.size,names=colnames(RGcpm),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Get log2 counts per million
logcounts <- cpm(RGcpm,log=TRUE)
# Check distributions of samples using boxplots
pdf('NonNormalizedLogCPM.pdf',w=30,h=10)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

##########################################################################
#full model with repeated measure correction 
#main effect of treatment: four groups
##################################################################################
#make design matrix
info$group <- as.character(info$group)
info$group <- gsub(" ","",info$group)
info$group2 <- paste(info$group,info$sex, sep="_")

info$group <- factor(info$group, levels=c("Saline_Brain","PolyIC_Brain", "Saline_Placenta", "PolyIC_Placenta"))

info$group2 <- factor(info$group2)

info$Dam <- factor(info$Dam)

info$sex <- as.factor(info$sex)

info$uterine.horn.position <- as.factor(info$uterine.horn.position)

#set design matrix
design <- model.matrix(~0 + group2 +uterine.horn.position, data=info) #if using a 0 intercept must set up contrasts

colnames(design)

#make contrasts
#https://support.bioconductor.org/p/91718/
#test whether the average treatment effect across all cell lines and treatment methods is significantly different from zero
#contrast compares the average of all treatment samples to the average of all control samples.

#https://seqqc.wordpress.com/2020/11/28/10-tips-tricks-for-complex-model-matrix-designs-in-dge-analysis/
#Contrasts for conditions:


# my.contrasts <- makeContrasts(
#   BrainSalvsPolyIC = groupPolyIC_Brain - groupSaline_Brain,
#   PlacentaSalvsPolyIC = groupPolyIC_Placenta - groupSaline_Placenta,
#   BrainSalvsPlacentaSal = groupSaline_Placenta - groupSaline_Brain,
#   BrainPolyICvsPlacentaPolyIC = groupPolyIC_Placenta - groupPolyIC_Brain,
#   levels=design)


my.contrasts <- makeContrasts(
  BrainSalvsPolyIC = (group2PolyIC_Brain_M + group2PolyIC_Brain_F)/2 - (group2Saline_Brain_M + group2Saline_Brain_F)/2,
  PlacentaSalvsPolyIC = (group2PolyIC_Placenta_M + group2PolyIC_Placenta_F)/2 - (group2Saline_Placenta_M + group2Saline_Placenta_F)/2,
  BrainSalvsPlacentaSal = (group2Saline_Placenta_M + group2Saline_Placenta_F)/2 - (group2Saline_Brain_M + group2Saline_Brain_F)/2,
  BrainPolyICvsPlacentaPolyIC = (group2PolyIC_Placenta_F + group2PolyIC_Placenta_M)/2 - (group2PolyIC_Brain_M + group2PolyIC_Brain_F)/2,
  
  M_BrainSalvsPolyIC = group2PolyIC_Brain_M - group2Saline_Brain_M,
  M_PlacentaSalvsPolyIC = group2PolyIC_Placenta_M - group2Saline_Placenta_M,
  M_BrainSalvsPlacentaSal = group2Saline_Placenta_M - group2Saline_Brain_M,
  M_BrainPolyICvsPlacentaPolyIC = group2PolyIC_Placenta_M - group2PolyIC_Brain_M,
  
  F_BrainSalvsPolyIC = group2PolyIC_Brain_F - group2Saline_Brain_F,
  F_PlacentaSalvsPolyIC = group2PolyIC_Placenta_F - group2Saline_Placenta_F,
  F_BrainSalvsPlacentaSal = group2Saline_Placenta_F - group2Saline_Brain_F,
  F_BrainPolyICvsPlacentaPolyIC = group2PolyIC_Placenta_F - group2PolyIC_Brain_F,
  
  levels=design)
#https://support.bioconductor.org/p/55070/



##################################################################################
#Normalization TMM and Voom
##################################################################################

#par(mfrow=c(1,2))
DGE=calcNormFactors(RGcpm,method =c("TMM")) #TMM normalization for library size/composition

pdf('VoomTrend.pdf',w=6,h=4)
v=voom(DGE,design,plot=T)
dev.off()

#account for repeated measures from same Dam
corfit <- duplicateCorrelation(v, design, block=info$Dam)
#corfit$consensus =  0.06360304

# # apply voom again (with the block and correlation parameters this time)
v <- voom(DGE, design, block = info$Dam, correlation = corfit$consensus)

fit <- lmFit(v, design, block = info$Dam, correlation = corfit$consensus)

fit2 <- contrasts.fit(fit,my.contrasts)
fit2 <- eBayes(fit2,robust=TRUE )



pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(fit2, main="Final model: Mean−variance trend")
dev.off()

#save
save(RG,RGcpm,info,design,my.contrasts,fit,fit2,v,file="DEGlist.RData")


#box plots for the voom normalised data to compare to before normalisation (only RLE)
#v$E already log2CPM
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="TMM & Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

#plot only normalized data
pdf('TMM_VOOM_NormalizedLogCPM.pdf',w=30,h=10)
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="TMM & Voom transformed logCPM")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
dev.off()



##################################################################################
#DE analysis:
##################################################################################
str(fit2)
#summary(decideTests(fit2,adjust.method="fdr", method="global"))
dt <- decideTests(fit2,adjust.method="BH", method="separate")
sum <- summary(dt)
#more weight to fold-changes in the ranking
#treat computes empirical Bayes moderated-t p-values relative to a minimum meaningful fold-change threshold. 
#Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true 
#log2-fold-change is greater than lfc in absolute value (McCarthy and Smyth, 2009). 
#In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc]. 
#When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and 
#for prioritizing genes that are biologically important. 
#treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods. 
#The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed. 
#When lfc=0, treat is identical to eBayes, except that F-statistics and B-statistics are not computed. 
#The lfc threshold is usually chosen relatively small, because significantly DE genes must all have fold changes substantially greater than the testing threshold. 
#Typical values for lfc are log2(1.1), log2(2) or log2(2). The top genes chosen by treat can be examined using topTreat.

# Treat relative to a ~2 fold-change
# tfit <- treat(fit2,lfc=log2(1))
# dt <- decideTests(tfit,adjust.method="BH", method="separate")
# sum <- summary(dt)

#no log fc filter:
write.csv(sum,"SummaryCount_DEGs.csv")

dt <- as.data.frame(dt)
tfit <- fit2



#######################################################################
#get out DE lists for each contrast:
#######################################################################
library("calibrate")
comparisons=(coef(tfit))
comparisons=colnames(comparisons)

comp_out <- as.data.frame(rownames(v$E))
names(comp_out) <- c("GeneID")

SumTableOut <- NULL

for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons 
  
  tmp=topTreat(tfit,coef=i,number=nrow(comp_out),adjust.method="BH")
  dim(tmp[(tmp$adj.P.Val<0.05),]) # number of DE genes
  
  #LogFC values:https://support.bioconductor.org/p/82478/
  tmp$direction <- c("none")
  tmp$direction[which(tmp$logFC > 0)] = c("Increase")
  tmp$direction[which(tmp$logFC < 0)] = c("Decrease")
  
  tmp$significance <- c("nonDE")
  tmp$significance[which(tmp$adj.P.Val <0.05)] = c("DE")
  
  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(tmp$significance,tmp$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)
  
  #get geneids  
  tmp$GeneID <- rownames(tmp)
  
  #gene gene names and expression levels
  tmp2 <- tmp
  
  tmp2$comparison <- paste(comp)
  
  #add in CPM values
  GeneSummaryCPM3 <- GeneSummaryCPM2[,c(1:5)]
  DECPM <- merge(tmp2, GeneSummaryCPM3,by.x="ensembl_gene_id", by.y="GeneID")
  
  write.csv(DECPM,file = paste(comp,"_DEgenes.csv"))
  
  #save to output:
  #merge <- merge(comp_out,tmp2, by= "GeneID")
  merge2 <- tmp2 %>% dplyr::select(GeneID,logFC,t,P.Value,adj.P.Val,direction,significance)
  colnames(merge2) <- paste(colnames(merge2),comp,sep=".")
  colnames(merge2)[1] <- c("GeneID")
  comp_out <- merge(comp_out, merge2, by="GeneID")
  
  #data for plot with gene names:
  genenames <- tmp2 %>% dplyr::select(adj.P.Val,logFC,mgi_symbol) %>% distinct()
  
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  #with(subset(genenames, logFC < -2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  #with(subset(genenames, logFC > 2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #color all significant genes
  with(subset(genenames, -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-2,2), col = "black", lty = 2, lwd = 1)
  
  # Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(genenames, adj.P.Val<0.05), textxy(logFC, -log10(adj.P.Val), labs=mgi_symbol, cex=.8))
  
  dev.off()
  
}


write.csv(SumTableOut,"SummaryTableDEgenes.csv")

#master outfile
mout <- merge(comp_out, GeneSummaryCPM2, by="GeneID")
write.csv(mout,"AllDEG_AllConditions.csv")

#######################################################################
#DE plots for Upregulated genes with LPS treatment
#######################################################################
library(VennDiagram)
#venn for LPS treatments:
#The number of genes that are not DE in either comparison are marked in the bottom-right.

pdf("Venn_PolyIC>Sal.pdf",height=10,width = 10)
vennDiagram(dt[,c("BrainSalvsPolyIC","PlacentaSalvsPolyIC")],
            circle.col=c("turquoise","purple"), include="up",
            names=c("Brain PolyIC>Sal","Placenta PolyIC>Sal"))
dev.off()


pdf("Venn_PolyIC<Sal.pdf",height=10,width = 10)
vennDiagram(dt[,c("BrainSalvsPolyIC","PlacentaSalvsPolyIC")],
            circle.col=c("turquoise","purple"), include="down",
            names=c("Brain PolyIC<Sal","Placenta PolyIC<Sal"))
dev.off()


pdf("Venn_allDEGs.pdf",height=10,width = 10)
vennDiagram(dt[,c("BrainSalvsPolyIC","PlacentaSalvsPolyIC")],
            circle.col=c("turquoise","purple"),
            names=c("Brain PolyICvsSal DEGs","Placenta PolyICvsSal DEGs"))
dev.off()



