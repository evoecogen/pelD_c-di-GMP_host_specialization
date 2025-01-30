#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("plotly")
#install.packages("VennDiagram")

library(VennDiagram)
library(ggplot2)
library(edgeR)
library(cowplot)
library(plotly)



#------------------------------------------------
# edgeR Workflow (R 4.4.1, EdgeR 4.2.1)
#------------------------------------------------
rm(list = ls())
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

my_theme <- function(data){
  
  # Set color palettes
  library(RColorBrewer)
  
  palette <- brewer.pal("Greys", n=9)
  col.back <-"grey75"
  col.grid <- palette[3]
  col.text <- palette[7]
  
  # Set plotting theme
  theme_bw(base_size=15) +
    
    # Set the entire chart region to a light gray color
    #theme(panel.background=element_rect(fill=col.back, color=col.back)) +
    #theme(plot.background=element_rect(fill=col.back, color=col.back)) +
    theme(panel.border=element_rect(fill= "transparent", color="black")) +
    theme(strip.background = element_rect(color="white",fill=col.back))+
    
    # Format the grid
    theme(panel.grid.major=element_line(color=col.grid,size=0.35)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Hide the legend
    #theme(legend.position="none") +
    theme(legend.background = element_rect(fill="transparent")) +
    theme(legend.text = element_text(size=12,color=col.text)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=col.text,size=10,vjust=1.25,hjust=0.5, face = "italic")) +
    theme(axis.text.x=element_text(size=14,color=col.text,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=14,color=col.text)) +
    theme(axis.title.x=element_text(size=14,color=col.text, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=col.text, vjust=1.25)) +
    theme(strip.text=element_blank())  }

# Import the raw mapped read data
counts <- read.delim("deseqFile")
colnames(counts) <- gsub(".rpkm_0","",colnames(counts))

# Import annotation and remove rRNA etc. from reads
counts$gene <- as.character(counts$gene)
annotation <- read.delim("Pseudomonas_lurida_MYb11_6243.ptt", skip = 2)
RNT <- read.delim("Pseudomonas_lurida_MYb11_6243.rnt", skip = 2)
all(counts[!(counts$gene %in% annotation$Synonym),]$gene == RNT$Synonym)
counts <- counts[counts$gene %in% annotation$Synonym,]
rownames(counts) <- counts$gene
counts$gene <- NULL
annotation <- annotation[,c(6,5,9,1,2,3)]
rownames(annotation) <- as.character(annotation$Synonym)

# Create grouping and sample info all
grouping_matrix <- matrix(0, nrow = 40, ncol = 4)
rownames(grouping_matrix) <- colnames(counts)
strings <- colnames(counts)
strings <- strsplit(strings, "_")
for(i in 1:length(colnames(counts))){
  grouping_matrix[i,] <- strings[[i]]
  grouping_matrix[i,4] <- paste(grouping_matrix[i,1], grouping_matrix[i,3], sep = ".")
}
grouping <- as.data.frame(grouping_matrix)
colnames(grouping) <- c("group","treatment", "replicate","strain")
grouping$group<-(as.factor(grouping$group))
grouping$strain<-(as.factor(grouping$strain))
grouping$treatment<-(as.factor(grouping$treatment))
grouping$replicate<-(as.factor(grouping$replicate))
write.csv(grouping, file="output/edgeR_norm/grouping.csv" )

# Create grouping and sample info liquid
groupingL <- droplevels(grouping[grouping$treatment == 'liquid',])
countsL<-counts
countsL<-countsL[grepl("_liquid", names(countsL))]
write.csv(groupingL, file="output/edgeR_norm/groupingL.csv" )

# Create grouping and sample info colony
groupingC <- droplevels(grouping[grouping$treatment == 'colony',])
countsC<-counts
countsC<-countsC[grepl("_colony", names(countsC))]
write.csv(groupingC, file="output/edgeR_norm/groupingC.csv" )

# Create DGEList object ALL
myList <- DGEList(counts = as.matrix(counts), group = grouping$group, genes = annotation)
keep <- filterByExpr(myList)
myList <- myList[keep, , keep.lib.sizes = FALSE]
myList <- calcNormFactors(myList)

# Create DGEList object liquid
myListL <- DGEList(counts = as.matrix(countsL), group = groupingL$group, genes = annotation)
keepL <- filterByExpr(myListL)
myListL <- myListL[keepL, , keep.lib.sizes = FALSE]
myListL <- calcNormFactors(myListL)

# Create DGEList object colony
myListC <- DGEList(counts = as.matrix(countsC), group = groupingC$group, genes = annotation)
keepC <- filterByExpr(myListC)
myListC <- myListC[keepC, , keep.lib.sizes = FALSE]
myListC <- calcNormFactors(myListC)

# PCA plot ALL
logCPM <- cpm(myList, log = TRUE)
pc = prcomp(t(logCPM))
ls(pc)
summary(pc)
#screeplot(pcL)
D_pca = as.data.frame(pc$x)
D_pca$group = paste(grouping$group)
D_pca$treatment = grouping$treatment
n = grouping$group
colors <- c("liquid" = "#5ab4ac", "colony"= "#d8b365")
shapes <- c("WT" =21, "wspF" = 22, "wspE" = 23, "rph" = 24)
p <- ggplot(D_pca,aes(x=PC1,y=PC2))
percentage <- round(pc$sdev^2 / sum(pc$sdev^2) * 100, 1); percentage <- paste( colnames(D_pca), "(", paste( as.character(percentage), "%", ")", sep="") )
pca1 <- p + stat_ellipse(aes(shape=group,color=treatment), lwd= 0.5)+geom_point(aes(shape=group, fill= treatment), size = 2, stroke= 0.25) + xlab(percentage[1]) + ylab(percentage[2]) + ggtitle("PCA") +
  coord_fixed()+
  labs(title = NULL , color = "Liquid", shape = "Strain") + 
  scale_shape_manual(values = shapes) + scale_fill_manual(values = colors) + scale_color_manual(values = colors)+
  my_theme()

pca2 <- plot_grid( pca1, ncol = 1)
pca2

# PCA plot liquid
logCPML <- cpm(myListL, log = TRUE)
pcL = prcomp(t(logCPML))
ls(pcL)
summary(pcL)
#screeplot(pcL)
D_pcaL = as.data.frame(pcL$x)
D_pcaL$group = paste(groupingL$group)
D_pcaL$treatment = groupingL$treatment
nL = groupingL$group
colors <- c("liquid" = "#5ab4ac")
shapes <- c("WT" =21, "wspF" = 22, "wspE" = 23, "rph" = 24)
pL <- ggplot(D_pcaL,aes(x=PC1,y=PC2))
percentageL <- round(pcL$sdev^2 / sum(pcL$sdev^2) * 100, 1); percentageL <- paste( colnames(D_pcaL), "(", paste( as.character(percentageL), "%", ")", sep="") )
pca1L <- pL + stat_ellipse(aes(color=group), size= 1)+geom_point(aes(shape=group, fill= treatment), size = 3, stroke= 1) + xlab(percentageL[1]) + ylab(percentageL[2]) + ggtitle("PCA") +
  coord_fixed()+
  labs(title = NULL , color = "Liquid", shape = "Strain") + 
  scale_shape_manual(values = shapes) + scale_fill_manual(values = colors) + scale_color_manual(values = c("#5ab4ac","#5ab4ac","#5ab4ac","#5ab4ac"))+
  theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "plain"), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "plain"), axis.text.y = element_text(size = 12))+
  theme_bw()

pca2L <- plot_grid( pca1L, ncol = 1, rel_heights = c(0.1, 1))
pca2L

# PCA plot colony
logCPMC <- cpm(myListC, log = TRUE)
pcC = prcomp(t(logCPMC))
ls(pcC)
summary(pcC)
#screeplot(pcC)
D_pcaC = as.data.frame(pcC$x)
D_pcaC$group = paste(groupingC$group)
D_pcaC$treatment = groupingC$treatment
nL = groupingC$group
colors <- c("colony" = "#d8b365")
shapes <- c("WT" =21, "wspF" = 22, "wspE" = 23, "rph" = 24)
pC <- ggplot(D_pcaC,aes(x=PC1,y=PC3))
percentageC <- round(pcC$sdev^2 / sum(pcC$sdev^2) * 100, 1); percentageC <- paste( colnames(D_pcaC), "(", paste( as.character(percentageC), "%", ")", sep="") )
pca1C <- pC + stat_ellipse(aes(color=group), size=1) + geom_point(aes(shape=group, fill= treatment), size = 3, stroke=1) + xlab(percentageC[1]) + ylab(percentageC[3]) + ggtitle("PCA") +
  labs(title = NULL , fill= "Treatment", shape = "Strain") + 
  scale_shape_manual(values = shapes) + scale_fill_manual(values = colors) + scale_color_manual(values = c("#d8b365","#d8b365","#d8b365","#d8b365")) +
  theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "plain"), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "plain"), axis.text.y = element_text(size = 12))+
  coord_fixed()+
  theme_bw()

pca2C <- plot_grid( pca1C, ncol = 1, rel_heights = c(0.1, 1) )
pca2C

#test cluster
distanceL <-dist(t(logCPML),method = "maximum")
clustersL <- hclust(distanceL, method ="average")
plot(clustersL)

distanceC <-dist(t(logCPMC),method = "maximum")
clustersC <- hclust(distanceC, method ="average")
plot(clustersC)

# Export logCPM and logRPKM
write.csv(logCPML, file = "output/edgeR_norm/Stat Tests/logCPML.csv")
lengthL <- myListL$genes$Length * 3
logRPKML <- rpkm(myListL, gene.length = lengthL, log = TRUE)
write.csv(logRPKML, file = "output/edgeR_norm/Stat Tests/logRPKML.csv")

write.csv(logCPMC, file = "output/edgeR_norm/Stat Tests/logCPMC.csv")
lengthC <- myListC$genes$Length * 3
logRPKMC <- rpkm(myListC, gene.length = lengthC, log = TRUE)
write.csv(logRPKMC, file = "output/edgeR_norm/Stat Tests/logRPKMC.csv")

# Design matrix and Dispersion estimation liquid
designL <- model.matrix(~0 + groupingL$group)
colnames(designL) <- levels(groupingL$group)
write.table(designL,file = "output/edgeR_norm/designL.txt")
myListL <- estimateDisp(myListL, design = designL, robust = TRUE)
plotBCV(myListL)
myListL$common.dispersion

# Design matrix and Dispersion estimation colony
designC <- model.matrix(~0 + groupingC$group)
colnames(designC) <- levels(groupingC$group)
write.table(designC,file = "output/edgeR_norm/designC.txt")
myListC <- estimateDisp(myListC, design = designC, robust = TRUE)
plotBCV(myListC)
myListC$common.dispersion

# Fitting Model liquid
fitL <- glmQLFit(myListL, design = designL, robust = TRUE)
head(fitL$coefficients)
plotQLDisp(fitL)

# Fitting Model colony
fitC <- glmQLFit(myListC, design = designC, robust = TRUE)
head(fitC$coefficients)
plotQLDisp(fitC)

# Testing for Differential Expression liquid
myContrastsL <- makeContrasts(LwspF = wspF - WT,
                             LwspE = wspE - WT,
                             Lrph = rph - WT,
                             levels=designL)

comparisonsL <- colnames(myContrastsL)

testL <- glmQLFTest(fitL, contrast = myContrastsL)
plotMD(testL)
abline(h=c(-1, 1), col="blue")

resultsL <- list()
for(i in 1:length(comparisonsL)){
  testL <- glmQLFTest(fitL, contrast = myContrastsL[,comparisonsL[i]])
  extractL <- topTags(testL, n = Inf, adjust.method = "BH")
  resultsL[[comparisonsL[i]]] <- extractL$table
}
calcL <- function(chrvec, fitL, myContrastsL){
  testL <- glmQLFTest(fitL, contrast = myContrastsL[,chrvec])
  extractL <- topTags(testL, n = Inf, adjust.method = "BH")
  return(extractL$table)
}
resultsL[["LwspF"]] <- calcL(chrvec = comparisonsL[1], fit = fitL, myContrasts = myContrastsL)
resultsL[["LwspE"]] <- calcL(chrvec = comparisonsL[2], fit = fitL, myContrasts = myContrastsL)
resultsL[["Lrph"]] <- calcL(chrvec = comparisonsL[3], fit = fitL, myContrasts = myContrastsL)

# Export CSVs
for(i in 1:length(names(resultsL))){
  write.csv(resultsL[[names(resultsL)[i]]], file = paste0("output/edgeR_norm/Stat Tests/", names(resultsL)[i], ".csv"))
}

# Testing for Differential Expression colony
myContrastsC <- makeContrasts(CwspF = wspF - WT,
                              CwspE = wspE - WT,
                              Crph = rph - WT,
                              levels=designC)

comparisonsC <- colnames(myContrastsC)

resultsC <- list()
for(i in 1:length(comparisonsC)){
  testC <- glmQLFTest(fitC, contrast = myContrastsC[,comparisonsC[i]])
  extractC <- topTags(testC, n = Inf, adjust.method = "BH")
  resultsC[[comparisonsC[i]]] <- extractC$table
}
calcC <- function(chrvec, fitC, myContrastsC){
  testC <- glmQLFTest(fitC, contrast = myContrastsC[,chrvec])
  extractC <- topTags(testC, n = Inf, adjust.method = "BH")
  return(extractC$table)
}
resultsC[["CwspF"]] <- calcC(chrvec = comparisonsC[1], fit = fitC, myContrasts = myContrastsC)
resultsC[["CwspE"]] <- calcC(chrvec = comparisonsC[2], fit = fitC, myContrasts = myContrastsC)
resultsC[["Crph"]] <- calcC(chrvec = comparisonsC[3], fit = fitC, myContrasts = myContrastsC)

# Export CSVs
for(i in 1:length(names(resultsC))){
  write.csv(resultsC[[names(resultsC)[i]]], file = paste0("output/edgeR_norm/Stat Tests/", names(resultsC)[i], ".csv"))
}

#------------------------------------------------
# Vulcanos and Venndiagramms
#------------------------------------------------

rm(list = ls())
my_theme <- function(data){
  
  # Set color palettes
  library(RColorBrewer)
  
  palette <- brewer.pal("Greys", n=9)
  col.back <-"grey75"
  col.grid <- palette[3]
  col.text <- palette[7]
  
  # Set plotting theme
  theme_bw(base_size=15) +
    
    # Set the entire chart region to a light gray color
    #theme(panel.background=element_rect(fill=col.back, color=col.back)) +
    #theme(plot.background=element_rect(fill=col.back, color=col.back)) +
    theme(panel.border=element_rect(fill= "transparent", color="black")) +
    theme(strip.background = element_rect(color="white",fill=col.back))+
    
    # Format the grid
    theme(panel.grid.major=element_line(color=col.grid,size=0.35)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Hide the legend
    #theme(legend.position="none") +
    theme(legend.background = element_rect(fill="transparent")) +
    theme(legend.text = element_text(size=12,color=col.text)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=col.text,size=10,vjust=1.25,hjust=0.5, face = "italic")) +
    theme(axis.text.x=element_text(size=14,color=col.text,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=14,color=col.text)) +
    theme(axis.title.x=element_text(size=14,color=col.text, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=col.text, vjust=1.25)) +
    theme(strip.text=element_blank())  }

#liquid
LwspF_liquid <- read.csv("output/edgeR_norm/Stat Tests/LwspF.csv")
LwspE_liquid <- read.csv("output/edgeR_norm/Stat Tests/LwspE.csv")
Lrph_liquid <- read.csv("output/edgeR_norm/Stat Tests/Lrph.csv")

arrange.table <- function(table){
  output <- table
  rownames(output) <- output$Synonym
  output$X <- NULL
  output <- output[order(output$Synonym),]
  return(output)
}

LwspF_liquid <- arrange.table(LwspF_liquid)
LwspE_liquid<- arrange.table(LwspE_liquid)
Lrph_liquid <- arrange.table(Lrph_liquid)


Lwrinkly <-merge(merge(LwspF_liquid, LwspE_liquid, by = c("Synonym","Gene","Product","Strand","Location","Length"),suffixes = c("wspF","wspE")), Lrph_liquid, by = c("Synonym","Gene","Product","Strand","Location","Length"))
write.csv(Lwrinkly, file = "output/edgeR_norm/Stat Tests/Lwrinkly_liquid.csv")

cans<-Lwrinkly[Lwrinkly$Synonym %in% c("CLM75_RS17855", "CLM75_RS17845","CLM75_RS01470","CLM75_RS02765","CLM75_RS19420", "CLM75_RS19415"),]
vulcano<-function(data){
  plot<-ggplot(data)+
  aes(y=-log10(FDRwspE), x=logFCwspE, text = paste("Symbol:",Synonym))+
  geom_point(aes(fill= FDRwspE<0.05 & (logFCwspE>1 | logFCwspE<(-1)), size= FDRwspE<0.05 & (logFCwspE>1 | logFCwspE<(-1)), colour = FDRwspE<0.05 & (logFCwspE>1 | logFCwspE<(-1))),shape = 21)+
    geom_point(data= cans, shape = 21, size = 2.5, fill = "#ef8a62", colour = "black") + 
    #geom_point(data = down_il_genes,
              # shape = 21,
               #size = 2, 
               #fill = "steelblue", 
               #colour = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey50", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#ef8a62", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#67a9cf", size=1) +
  geom_label(data = cans, aes(label = Synonym), nudge_y = -0.25) +
  scale_fill_manual(values = c("grey90","grey40")) +
  scale_colour_manual(values = c("grey90","black")) +
  scale_size_manual(values = c(1,2)) +
  annotate("rect", xmin = -3, xmax = 5, ymin = 1, ymax = 10, alpha=0) +
    labs(title="Volcano plot", caption=paste0("produced on ", Sys.time())) +
    my_theme()+
  theme(legend.title = element_blank())
return(plot)}

#Vulcanoplots
vulcano(Lwrinkly)

#colony
CwspF_colony <- read.csv("output/edgeR_norm/Stat Tests/CwspF.csv")
CwspE_colony <- read.csv("output/edgeR_norm/Stat Tests/CwspE.csv")
Crph_colony <- read.csv("output/edgeR_norm/Stat Tests/Crph.csv")

CwspF_colony <- arrange.table(CwspF_colony)
CwspE_colony<- arrange.table(CwspE_colony)
Crph_colony <- arrange.table(Crph_colony)

Cwrinkly <-merge(merge(CwspF_colony, CwspE_colony, by = c("Synonym","Gene","Product","Strand","Location","Length"),suffixes = c("wspF","wspE")), Crph_colony, by = c("Synonym","Gene","Product","Strand","Location","Length"))
write.csv(Cwrinkly, file = "output/edgeR_norm/Stat Tests/Cwrinkly_colony.csv")

cans<-Cwrinkly[Cwrinkly$Synonym %in% c("CLM75_RS05020", "CLM75_RS05035"),]
vulcano(Cwrinkly)

#--------------------------------------------
# Introduction of Cutoff
#--------------------------------------------
rm(list = ls())

#liquid
LwspF_liquid <- read.csv("output/edgeR_norm/Stat Tests/LwspF.csv")
LwspE_liquid <- read.csv("output/edgeR_norm/Stat Tests/LwspE.csv")
Lrph_liquid <- read.csv("output/edgeR_norm/Stat Tests/Lrph.csv")

arrange.table <- function(table){
  output <- table
  rownames(output) <- output$Synonym
  output$X <- NULL
  output <- output[order(output$Synonym),]
  return(output)
}
LwspF_liquid <- arrange.table(LwspF_liquid)
LwspE_liquid <- arrange.table(LwspE_liquid)
Lrph_liquid <- arrange.table(Lrph_liquid)

select.sig <- function(table, cutfdr, cutfc){
  anydiff <- (abs(table[,7]) >= cutfc)
  tablered <- table[anydiff,]
  tablered <- tablered[tablered$FDR < cutfdr,]
  #tablered$FDR <- p.adjust(tablered$PValue, method = "BH")
  tablered <- tablered[order(tablered$FDR, tablered$PValue),]
  return(tablered)
}
LwspF_liquid_fdr <- select.sig(table = LwspF_liquid, cutfdr = 0.05, cutfc = 0)
LwspE_liquid_fdr <- select.sig(table = LwspE_liquid, cutfdr = 0.05, cutfc = 0)
Lrph_liquid_fdr <- select.sig(table = Lrph_liquid, cutfdr = 0.05, cutfc = 0)

write.csv(LwspF_liquid_fdr, file = "output/edgeR_norm/Stat Tests/LwspF_liquid_fdr.csv")
write.csv(LwspE_liquid_fdr, file = "output/edgeR_norm/Stat Tests/LwspE_liquid_fdr.csv")
write.csv(Lrph_liquid_fdr, file = "output/edgeR_norm/Stat Tests/Lrph_liquid_fdr.csv")

#Venndiagramm
venndata<-list(rph=Lrph_liquid_fdr$Synonym,wspF=LwspF_liquid_fdr$Synonym,wspE=LwspE_liquid_fdr$Synonym)
venn.plot<-venn.diagram(venndata,filename=NULL,main="rph vs. wspF vs. wspE",compression = "lzw",lwd = 2 )
grid.newpage()
grid.draw(venn.plot)

Lwrinkly_sig_fdr <-merge(merge(LwspF_liquid_fdr, LwspE_liquid_fdr, by = c("Synonym","Gene","Product","Strand","Location","Length"),suffixes = c("wspF","wspE")), Lrph_liquid_fdr, by = c("Synonym","Gene","Product","Strand","Location","Length"))
write.csv(Lwrinkly_sig_fdr, file = "output/edgeR_norm/Stat Tests/Lwrinkly_liquid_sig_fdr.csv")

#colony
CwspF_colony <- read.csv("output/edgeR_norm/Stat Tests/CwspF.csv")
CwspE_colony <- read.csv("output/edgeR_norm/Stat Tests/CwspE.csv")
Crph_colony <- read.csv("output/edgeR_norm/Stat Tests/Crph.csv")

CwspF_colony <- arrange.table(CwspF_colony)
CwspE_colony <- arrange.table(CwspE_colony)
Crph_colony <- arrange.table(Crph_colony)

CwspF_colony_fdr <- select.sig(table = CwspF_colony, cutfdr = 0.05, cutfc = 0)
CwspE_colony_fdr <- select.sig(table = CwspE_colony, cutfdr = 0.05, cutfc = 0)
Crph_colony_fdr <- select.sig(table = Crph_colony, cutfdr = 0.05, cutfc = 0)

write.csv(CwspF_colony_fdr, file = "output/edgeR_norm/Stat Tests/CwspF_colony_fdr.csv")
write.csv(CwspE_colony_fdr, file = "output/edgeR_norm/Stat Tests/CwspE_colony_fdr.csv")
write.csv(Crph_colony_fdr, file = "output/edgeR_norm/Stat Tests/Crph_colony_fdr.csv")

#Venndiagram
Cwrinkly_colony_sig_fdr <-merge(merge(CwspF_colony_fdr, CwspE_colony_fdr, by = c("Synonym","Gene","Product","Strand","Location","Length"),suffixes = c('wspF','wspE')), Crph_colony_fdr, by = c("Synonym","Gene","Product","Strand","Location","Length"),suffixes='rph')

venndata<-list(rph=Crph_colony_fdr$Synonym,wspF=CwspF_colony_fdr$Synonym,wspE=CwspE_colony_fdr$Synonym)
venn.plot<-venn.diagram(venndata,filename=NULL,main="rph vs. wspF vs. wspE")
grid.newpage()
grid.draw(venn.plot)

write.csv(Cwrinkly_colony_sig_fdr, file = "output/edgeR_norm/Stat Tests/Cwrinkly_colony_sig_fdr.csv")