#----------------------------------
# 2024
# in vitro biofilm
#----------------------------------

#install.packages('openxlsx')
#install.packages("rstatix")

library(ggplot2)
library(openxlsx)
library(reshape2)
library(rstatix)
library(lme4)
library(multcomp)

# Setting the directory for reading in data
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

#--- (1) Collect data -------------------------------

data1<-read.xlsx("Crystal_violet_01082024.xlsx")
data2<-read.xlsx("Crystal_violet_07082024.xlsx")
data3<-read.xlsx("Crystal_violet_15082024_3.xlsx")
data4<-read.xlsx("Crystal_violet_15082024_4.xlsx")

# Subset data of interest
#Absorbance 550 nm
SubdataA<-function(data)
{
dataA550 <- data[55:63, 1:13]

# Create 96-well-plate layout
dataA <- as.matrix(dataA550[2:9,2:13])
rownames(dataA) <- dataA550[2:9,1]
colnames(dataA) <- dataA550[1,2:13]

# Transform matrix into column format
dataA<- melt(dataA)
well<-paste(dataA$Var1,dataA$Var2,sep='')
tmp<-data.frame(well,value=dataA$value);tmp

dataRep <-tmp
return(dataRep)
}

SubdataA2<-function(data)
{
  dataA550 <- data[19:27, 1:13]
  
  # Create 96-well-plate layout
  dataA <- as.matrix(dataA550[2:9,2:13])
  rownames(dataA) <- dataA550[2:9,1]
  colnames(dataA) <- dataA550[1,2:13]
  
  # Transform matrix into column format
  dataA<- melt(dataA)
  well<-paste(dataA$Var1,dataA$Var2,sep='')
  tmp<-data.frame(well,value=dataA$value);tmp
  
  dataRep <-tmp
  return(dataRep)
}

dataRep1<- SubdataA(data1)
dataRep2<- SubdataA2(data2)
dataRep3<- SubdataA(data3)
dataRep4<- SubdataA(data4)

# Add data legend
addLegend<-function(data)
{
data <- merge(data, legend, by = "well")
data <- merge(data, legend2, by = "id")
data$id <- as.factor(data$id)
data$bacteria <- as.factor(data$bacteria)
data$value <- as.numeric(data$value)
return(data)
}

legend <- read.csv("layout_1.csv")
legend2 <- read.csv("strainIDs.csv", header = T,sep = ";"  )
dataRep1<-addLegend(dataRep1)
rm(legend)

legend <- read.csv("layout_2.csv",sep = ";"  )
dataRep2<-addLegend(dataRep2)
rm(legend)

legend <- read.csv("layout_3.csv",sep = "," )
dataRep3<-addLegend(dataRep3)
rm(legend)

legend <- read.csv("layout_4.csv",sep = ";")
dataRep4<-addLegend(dataRep4)
rm(legend)

#--- (2) Determine absorbance 550nm data -------------------------------

Ab<-function(data)
{
  
#mean absorcance 550nm
#blank raw data
blank <- data[data$bacteria == 'TSB',]
blank<-subset(blank, select = value,)
blank$mean<-mean(blank$value)
blank<-blank[!duplicated(blank$mean),]
data$value<-data$value-blank$mean

#mean absorbance per technical replicate
data$value<-ave(data$value, data$bacteria, data$id, FUN = mean)
data<-subset(data, select = -well)
data<-data[!duplicated(data),]

return(data)
}

dataRep1 <- Ab(dataRep1)
dataRep1$rep <-1
dataRep1$rep<-as.factor(dataRep1$rep)

dataRep2 <- Ab(dataRep2)
dataRep2$rep <-2
dataRep2$rep<-as.factor(dataRep2$rep)

dataRep3 <- Ab(dataRep3)
dataRep3$rep <-3
dataRep3$rep<-as.factor(dataRep3$rep)

dataRep4 <- Ab(dataRep4)
dataRep4$rep <-4
dataRep4$rep<-as.factor(dataRep4$rep)

dataRep <-rbind(dataRep1,dataRep2,dataRep3,dataRep4)

#remove TSB 
dataRep<-dataRep[dataRep$bacteria !='TSB',]

#--- (3) Plot data -------------------------------
####9. Boxplot ####

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
    theme(legend.position="none") +
    #theme(legend.background = element_rect(fill=color.background)) +
    #theme(legend.text = element_text(size=12,color=col.text)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=col.text,size=10,vjust=1.25,hjust=0.5, face = "italic")) +
    theme(axis.text.x=element_text(size=14,color=col.text,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=14,color=col.text)) +
    theme(axis.title.x=element_text(size=14,color=col.text, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=col.text, vjust=1.25)) +
    theme(strip.text=element_blank())  }

# ABsorbance 550nm
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}

PlotAb<- function(Data){
  #Data$strain <- factor(Data$strain,levels = c('WT','MT12','MT21','MT25','MT26','MT13','MT33','MT11'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=bacteria,y=value))+
    geom_boxplot(aes(),lwd=0.5,fatten=1, colour= "black",outlier.shape = NA)+
    geom_point(aes(linetype = "Replicate"),shape=21,size=2,fill="grey90", colour="black", position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0,seed=NA))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y="Absorbance 550 nm")
  
  return(Pop_Plot)
}

dataRep$bacteria <- factor(dataRep$bacteria,levels = c('OP50','MYb11','MT14','alg','alg44','azu','fimA','fli','usher','pel','fimA usher','fli usher','pel usher'),ordered = TRUE)
PlotAb (dataRep)

#define comparisons of interest
comp48<-cbind.data.frame(rep("MYb11", 12), c('MT14','alg','alg44','azu','fimA','fli','usher','pel','pel usher','fimA usher','fli usher', 'OP50'))
comp14<-cbind.data.frame(rep("MT14", 12), c('MYb11','alg','alg44','azu','fimA','fli','usher','pel','pel usher','fimA usher','fli usher', 'OP50'))

stat<-function(comparisons,data,test){
  
  # Create object with test results
  res           <- as.data.frame(matrix(nrow = nrow(comparisons), ncol = 17))
  colnames(res) <- c("Group 1", "Group 2", "Shapiro-Wilk test statistic Group 1", "Shapiro-Wilk p value Group 1", "Shapiro-wilk test statistic Group 2", "Shapiro-Wilk p value Group 2", "F test: F statistic", "F test: num df", "F test: denom df", "F test: pvalue", "test for comparisons of central tendencies", "H0", "test statistic", "df", "p value", "FDR corrected p value", "reject H0")
  
  for (j in 1:nrow(res)) {
    # extract readings for group 1 and group 2
    data1       <- data$value[data$bacteria == comparisons[j, 1] ]
    data2       <- data$value[data$bacteria == comparisons[j, 2] ]
    
    # perform shapiro tests
    shapiro1    <- shapiro.test(data1)
    shapiro2    <- shapiro.test(data2)
    
    # perform F test
    f.test      <- var.test(data1, data2)
    
    res[j, 1]   <- comparisons[j, 1]
    res[j, 2]   <- comparisons[j, 2]
    res[j, 3]   <- shapiro1$statistic
    res[j, 4]   <- shapiro1$p.value
    res[j, 5]   <- shapiro2$statistic
    res[j, 6]   <- shapiro2$p.value
    res[j, 7]   <- f.test$statistic
    res[j, 8:9] <- f.test$parameter
    res[j, 10]  <- f.test$p.value
    res[j, 11]  <- ifelse(shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05 & f.test$p.value > 0.05, "t test with equal variances", ifelse(shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05 & f.test$p.value <= 0.05, "Welch's t test", "Mann-Whitney U test"))
    res[j, 12]  <- ifelse(res[j, 11] == "Mann-Whitney U test", "median group 1 = medan group 2", "mean group 1 = mean group 2")
    
    # perform test for comparisons of central tendencies
    cent.test   <- NULL
    if (res[j, 11] == "t test with equal variances") {
      cent.test <- t.test(data1, data2, var.equal = T)
    }
    else if (res[j, 11] == "Welch's t test") {
      cent.test <- t.test(data1, data2, var.equal = F)
    }
    else {
      cent.test <- wilcox.test(data1, data2)
    }
    
    res[j, 13]  <- cent.test$statistic
    res[j, 14]  <- ifelse(is.null(cent.test$parameter), NA, cent.test$parameter)
    res[j, 15]  <- cent.test$p.value
  }
  
  # perform FDR correction
  res[, 16]    <- p.adjust(res[, 15], method = "fdr")
  res[, 17]    <- ifelse(res[, 16] < 0.05, "yes", "no")
  
  
  # save results as csv-file
  write.csv(res,file = paste('Stats_out/stats_',test,'.csv',sep=''),quote = F,row.names = F)
}

#execute function
stat(comp48,dataRep,"anc")
stat(comp14,dataRep,"wspE")

stats_anc<-read.csv("Stats_out/stats_anc.csv",TRUE,",")
stats_wspE<-read.csv("Stats_out/stats_wspE.csv",TRUE,",")