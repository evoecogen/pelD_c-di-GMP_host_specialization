#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#install.packages("ggplot2")
#install.packages("rstudioapi")
#install.packages("tidyverse")

library(clusterProfiler) #version 4.12.6
library(ggplot2)
library(dplyr)

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


#-----------------------------------------
# Data Import and Database construction
#-----------------------------------------
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)

wspF_c <- read.csv2("output/edgeR_norm/Stat Tests/CwspF.csv")
rownames(wspF_c) <- wspF_c$X
wspF_c$X <- NULL


# ANOVA Results

wrinkly_cf <- read.csv2("output/edgeR_norm/Stat Tests/Cwrinkly_colony_sig_fdr.csv")
rownames(wrinkly_cf) <- wrinkly_cf$X
wrinkly_cf$X <- NULL

wrinkly_lf <- read.csv2("output/edgeR_norm/Stat Tests/Lwrinkly_liquid_sig_fdr.csv")
rownames(wrinkly_lf) <- wrinkly_lf$X
wrinkly_lf$X <- NULL

# Define Background based on old locus tags
background <- as.character(wspF_c$Synonym)
background <- sort(background)

# Gene description data frame
descr <- data.frame(locus = wspF_c$Synonym, description = wspF_c$Product)
descr <- descr[order(descr$locus),]

# Data provided by pseudomonas.com
# Define GO Data
GO <- read.csv("gene_ontology_csv.csv")
GO.reduced <- GO[,c(1,2,3,5,6,7)]
GO.reduced <- unique(GO.reduced)

GO.MF <- GO.reduced[GO.reduced$Namespace == "molecular_function",]
GO.BP <- GO.reduced[GO.reduced$Namespace == "biological_process",]
GO.CC <- GO.reduced[GO.reduced$Namespace == "cellular_component",]

GO.MF$Accession<-as.factor(GO.MF$Accession)
GO.MF$GO.Term<-as.factor(GO.MF$GO.Term)
GO.MF$Locus.Tag<-as.factor(GO.MF$Locus.Tag)

GO.BP$Accession<-as.factor(GO.BP$Accession)
GO.BP$GO.Term<-as.factor(GO.BP$GO.Term)
GO.BP$Locus.Tag<-as.factor(GO.BP$Locus.Tag)

GO.CC$Accession<-as.factor(GO.CC$Accession)
GO.CC$GO.Term<-as.factor(GO.CC$GO.Term)
GO.CC$Locus.Tag<-as.factor(GO.CC$Locus.Tag)

GO.reduced$Accession<-as.factor(GO.reduced$Accession)
GO.reduced$GO.Term<-as.factor(GO.reduced$GO.Term)
GO.reduced$Locus.Tag<-as.factor(GO.reduced$Locus.Tag)

make.2gene <- function(GO.component){
  c2gene <- GO.component[,c(4,1)]
  c2gene <- c2gene[order(c2gene$Accession, c2gene$Locus.Tag),]
  c2gene$Accession <- droplevels(c2gene$Accession)
  c2gene$Locus.Tag <- droplevels(c2gene$Locus.Tag)
  rownames(c2gene) <- NULL
  return(c2gene)
}
make.2name <- function(GO.component){
  c2name <- GO.component[,4:5]
  c2name <- unique(c2name)
  c2name <- c2name[order(c2name$Accession),]
  c2name$Accession <- droplevels(c2name$Accession)
  c2name$GO.Term <- droplevels(c2name$GO.Term)
  rownames(c2name) <- NULL
  return(c2name)
}

MF2gene <- make.2gene(GO.MF)
MF2name <- make.2name(GO.MF)

BP2gene <- make.2gene(GO.BP)
BP2name <- make.2name(GO.BP)

CC2gene <- make.2gene(GO.CC)
CC2name <- make.2name(GO.CC)

GO2gene <- make.2gene(GO.reduced)
GO2name <- make.2name(GO.reduced)

#--------------------------------------------
# Gene Set Enrichment Analysis (GSEA)
#--------------------------------------------

#mix
make.geneList <- function(comp){
  fcs <- comp[,c(7,12,17)]
  rownames(fcs) <- comp[,1]
  means <- c()
  for(i in 1:nrow(fcs)){
    mean <- mean(c(fcs[i,1], fcs[i,2], fcs[i,3]))
    means <- c(means, mean)
  }
  names(means) <- rownames(fcs)
  means <- sort(means, decreasing = TRUE)
  return(means)
}

wrinkly_lf.list <- make.geneList(wrinkly_lf)
wrinkly_cf.list <- make.geneList(wrinkly_cf)

set.readable <- function(enrich.df, descr){
  column <- character()
  for(i in 1:nrow(enrich.df)){
    split <- strsplit(enrich.df$core_enrichment[i], split = "/")
    unlist <- unlist(split)
    descriptions <- descr$description[match(unlist, descr$locus)]
    merged <- paste(descriptions, collapse = "/")
    column <- c(column, merged)
  }
  return(column)
}

geneLists <- list(wrinkly_cf = wrinkly_cf.list,
                  wrinkly_lf = wrinkly_lf.list)
# Gene Ontology
namespace <- names(geneLists)
enrich.results <- list()
for(i in 1:length(namespace)){
  geneList <- geneLists[[namespace[i]]]
  set.seed(20240924)
  MFenrich <- GSEA(geneList = geneList, nPermSimple = 1000, minGSSize = 2,  pvalueCutoff = 1, pAdjustMethod = "fdr",
                   TERM2GENE = MF2gene, TERM2NAME = MF2name)
  MFenrich <- as.data.frame(MFenrich)
  MFenrich <- cbind("category" = rep("MF", nrow(MFenrich)), MFenrich)
  set.seed(20240924)
  BPenrich <-  GSEA(geneList = geneList, nPermSimple = 1000, minGSSize = 2, pvalueCutoff = 1, pAdjustMethod = "fdr",
                   TERM2GENE = BP2gene, TERM2NAME = BP2name)
  BPenrich <- as.data.frame(BPenrich)
  BPenrich <- cbind("category" = rep("BP", nrow(BPenrich)), BPenrich)
  set.seed(20240924)
  CCenrich <- GSEA(geneList = geneList, nPermSimple = 1000, minGSSize = 2,  pvalueCutoff = 1, pAdjustMethod = "fdr",
                   TERM2GENE = CC2gene, TERM2NAME = CC2name)
  CCenrich <- as.data.frame(CCenrich)
  CCenrich <- cbind("category" = rep("CC", nrow(CCenrich)), CCenrich)
  
  enrich.df <- rbind(MFenrich, BPenrich, CCenrich)
  if(nrow(enrich.df) >=1){enrich.df$gene.descriptions <- set.readable(enrich.df = enrich.df, descr = descr)}
  enrich.results[[namespace[i]]] <- enrich.df
}
# Export CSVs
for(i in 1:length(enrich.results)){
  write.csv2(enrich.results[[i]], file = paste("output/edgeR_norm/GSEA/GO_", names(enrich.results)[i], ".csv", sep = ""))
}

#plot GSEA test

GSEA_wcf <- read.csv2("output/edgeR_norm/GSEA/GO_wrinkly_cf.csv")
GSEA_wlf <- read.csv2("output/edgeR_norm/GSEA/GO_wrinkly_lf.csv")

# add a variable to this result that matches enrichment direction with phenotype
GSEA_plot<-function(data,title){
  
data <- data %>%
  mutate(Transcription = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"))

# create 'bubble plot' to summarize y signatures across x phenotypes
data$Transcription<-as.factor(data$Transcription)
data$Transcription<-factor(data$Transcription,levels = c("Up","Down"))
  data$core_enrichment<-as.factor(data$core_enrichment)
  data<-data[!duplicated(data$core_enrichment),]
  data<- na.omit(data)
  data<-droplevels(data)
  levels(data$Description)<- c("calcium ion binding","transmembrane transport","catalytic activity","phosphopantetheine binding", 
   "regulation of transcription, DNA-templated",
   "cell wall organization/chaperone-mediated protein folding/pilus organization/outer membrane-bounded periplasmic space",
   "cell adhesion/pilus")
data$Description <-reorder(data$Description, data$NES, increasing = T)
plot<-ggplot(data, aes(x=Transcription, y=Description)) + 
  geom_point(aes(size= setSize, color = p.adjust)) +
  scale_color_gradient(high= "#67a9cf", low = "#ef8a62")+
  my_theme()+
  theme(aspect.ratio=7/1)+
  labs(title= title,
       caption=paste0("produced on ", Sys.time()))
return(plot)}

GSEA_plot(GSEA_wlf[GSEA_wlf$p.adjust<=0.05,], title = "LW")

GSEA_plot(GSEA_wcf[GSEA_wcf$p.adjust<=0.05,], title= "CW")

