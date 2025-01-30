#########################################
# Short-term persistence competition    # 
######################################### 

#### get packages ####

#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("DescTools")

#### load packages####

library(rstudioapi)
library(ggplot2)
library(multcomp)

#### get name of the directory this file is in####

currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)

#### set working directory to currentDirectory####

setwd(currentDirectory)

#### get rawdata ####

Competition <-read.csv("202404_STP.csv",TRUE,";")
Competition2 <-read.csv("202407_STP.csv",TRUE,";")

#### function for rawdata ####

Rawdata <- function(Data){
# rename colume
names(Data)
names(Data)[1]<-"Replicate"
names(Data)[5]<-"1"
names(Data)[6]<-"2"
names(Data)[7]<-"3"
# edit data
Data$MT<-as.factor(Data$MT)
Data$Replicate<-as.factor(Data$Replicate)
Data$Set<-as.factor(Data$Set)
return(Data)
}

#### apply function####

Competition<-Rawdata(Competition)
Competition2<-Rawdata(Competition2)

#### function for CFUs/Worm ####

CFUSCompCol<-function(Data){
  
  # Subset two data frames, one with only worm data and one with only supernatants
  onlyWorms <- Data[Data$Typ == "W", ]
  onlySup <- Data[Data$Typ == "S", ]
  
  # Merge worm and sup data frames
  combiWormSup <- merge(onlyWorms, onlySup, by = c("Replicate", "MT","Set","Worms","Sample","Col"),
                        all.x = T)
  
  # Mean
  combiWormSup$MeanW <- rowMeans(combiWormSup[(8:10)],TRUE)
  combiWormSup$MeanS <- rowMeans(combiWormSup[(13:15)],TRUE)
  
  # Mean/Worm
  combiWormSup$'MeanW/Worm' <- (combiWormSup$'MeanW'/combiWormSup$'Worms')
  combiWormSup$'MeanS/Worm' <- (combiWormSup$'MeanS'/combiWormSup$'Worms')
  
  # CFU/5?l/Worm
  combiWormSup$'MeanW/Worm/25?l' <- (combiWormSup$'MeanW/Worm'* 10^combiWormSup$'Dilution.x')
  combiWormSup$'MeanS/Worm/25?l' <- (combiWormSup$'MeanS/Worm'* 10^combiWormSup$'Dilution.y')
  
  # CFU/100?l/Worm
  combiWormSup$'MeanW/Worm/100?l' <- (combiWormSup$'MeanW/Worm/25?l'*4)
  combiWormSup$'MeanS/Worm/100?l' <- (combiWormSup$'MeanS/Worm/25?l'*4)
 
  # CFUs
  combiWormSup$'CFUs' <- (combiWormSup$'MeanW/Worm/100?l'-combiWormSup$'MeanS/Worm/100?l')
 
  #final Tabel
  CFUs <- combiWormSup[,c(1:5,25)]
  CFUs2<-log10(CFUs[6])
  CFUs<-cbind(CFUs,CFUs2)
  names(CFUs)[7]<-"LOG_CFU"
  return(CFUs)
}

#### apply function####

CFUs_Comp <- CFUSCompCol(Competition)
CFUs_Comp2 <- CFUSCompCol(Competition2)

#### Prepare data for Box-Plot Fig. S4, 2a####

#CFUs_Comp$CFUs[CFUs_Comp$CFUs == 0]<-0.001
CFUs_Comp$CFUs[CFUs_Comp$CFUs< 0]<-NA
CFUs_Comp<- CFUs_Comp[!is.na(CFUs_Comp$CFUs),]
CFUs_Comp$MT<-as.factor(CFUs_Comp$MT)

CFUs_Comp2$CFUs[CFUs_Comp2$CFUs< 0]<-NA
CFUs_Comp2<- CFUs_Comp2[!is.na(CFUs_Comp2$CFUs),]
CFUs_Comp2$MT<-as.factor(CFUs_Comp2$MT)

#### Boxplot ####

scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}

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

CFUs_Comp<- CFUs_Comp[!is.infinite(CFUs_Comp$LOG_CFU),]
CFUs_Comp2<- CFUs_Comp2[!is.infinite(CFUs_Comp2$LOG_CFU),]

CFUs_Comp$Round<-as.factor(1)
CFUs_Comp2$Round<-as.factor(2)
STP_Solo<-rbind(CFUs_Comp,CFUs_Comp2)
STP_Solo$Set <- factor(STP_Solo$Set,levels = c('dT','MT48','MT14','algD','alg44','azu','fimA','fli', 'usher', 'pel', 'fimA/usher', 'fli/usher', 'pel/usher'),ordered = TRUE)
STP_Solo <-na.omit(STP_Solo)

#Fig. S4
PlotSTP_Solo<- function(Data){
  Pop_Plot <- ggplot(Data,aes(x=Set,y=CFUs))+
    geom_boxplot(lwd=0.5,fatten=1, colour= "black",outlier.shape = NA)+
    #theme(legend.title=element_text(size=12),legend.text=element_text(size=11))+
    geom_point(aes(linetype = "Replicate"),shape=21,size=2,fill="grey90", colour="black", position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0,seed=NA))+
    scale_y_log10(limits=c(1, 100000), breaks = c( 10, 100, 1000, 10000, 100000), labels = scientific_10 )+
    facet_grid(cols = vars(Round), scales = "free",space = "free")+
    my_theme()+
    labs(y="Short-term persistence [Mutant/MYb11dT]")
  
  return(Pop_Plot)
}

#Fig. S4
PlotSTP_Solo(STP_Solo[STP_Solo$Set %in% c("MT14","dT","azu","alg44", "algD", "fimA", "fli", "usher", "MT48", "pel", "pel/usher", "fimA/usher", "fli/usher"),])

#Fig. 2a
#subset data relative competition
Pup<- CFUs_Comp
MTs <- Pup[!Pup$MT == "dT", ]
dT<-Pup[Pup$MT== "dT",]
Pup <- merge(MTs, dT, by = c("Replicate","Set","Worms","Sample"),all.x = T)

Pup$'Rel' <- (Pup$'CFUs.x'/Pup$'CFUs.y')
Pup <- Pup[!is.na(Pup$MT.y),]

Pup$Set <- as.factor(Pup$Set)

Pup2<- CFUs_Comp2
MTs2 <- Pup2[!Pup2$MT == "dT", ]
dT2<-Pup2[Pup2$MT== "dT",]
Pup2 <- merge(MTs2, dT2, by = c("Replicate","Set","Worms","Sample"),all.x = T)

Pup2$'Rel' <- (Pup2$'CFUs.x'/Pup2$'CFUs.y')
Pup2 <- Pup2[!is.na(Pup2$MT.y),]

Pup2$Set <- as.factor(Pup2$Set)

Pup$Set <- factor(Pup$Set,levels = c('MT48+dT','MT14+dT','algD+dT','alg44+dT','azu+dT','fimA+dT','fli+dT','usher+dT'),ordered = TRUE)
Pup2$Set <- factor(Pup2$Set,levels = c('MT48+dT','MT14+dT','pel+dT','fimA/usher+dT','fli/usher+dT','pel/usher+dT'),ordered = TRUE)

Pup$Round<-as.factor(1)
Pup2$Round<-as.factor(2)
STP<-rbind(Pup,Pup2)

#Fig. 2a
PlotSTP<- function(Data){
  Pop_Plot <- ggplot(Data,aes(x=Set,y=Rel))+
    geom_hline(yintercept = 1,linetype = "dashed",size=0.5)+
    geom_boxplot(lwd=0.5,fatten=1, colour= "black",outlier.shape = NA)+
    #theme(legend.title=element_text(size=12),legend.text=element_text(size=11))+
    geom_point(aes(linetype = "Replicate"),shape=21,size=2,fill="grey90", colour="black", position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0,seed=NA))+
    scale_y_log10(limits=c(0.1, 100),labels = scientific_10)+
    facet_grid(cols = vars(Round), scales = "free",space = "free")+
    my_theme()+
    #theme(text = element_text(size=12),strip.text.x= element_text(colour="black"),strip.background = element_rect(colour="black", fill="white"),panel.grid = element_line(color = "grey", linetype = "dashed"),axis.title.x = element_blank(),axis.text.x=element_text(angle=15,hjust=1))+
    labs(y="Short-term persistence [Mutant/MYb11dT]")
  
  return(Pop_Plot)
}

#Fig. 2a
PlotSTP(STP)

###Total CFUs###
#### function for CFUs/Worm Fig. 2b####
#Fig. 2b

CFUSCompColT<-function(Data){
  
  # Subset two data frames, one with only worm data and one with only supernatants
  onlyWorms <- Data[Data$Typ == "W", ]
  onlySup <- Data[Data$Typ == "S", ]
  
  # Merge worm and sup data frames
  combiWormSup <- merge(onlyWorms, onlySup, by = c("Replicate", "MT","Set","Worms","Sample","Col"),
                        all.x = T)
  
  # Mean
  combiWormSup$MeanW <- rowMeans(combiWormSup[(8:10)],TRUE)
  combiWormSup$MeanS <- rowMeans(combiWormSup[(13:15)],TRUE)
  
  # Mean/Worm
  combiWormSup$'MeanW/Worm' <- (combiWormSup$'MeanW'/combiWormSup$'Worms')
  combiWormSup$'MeanS/Worm' <- (combiWormSup$'MeanS'/combiWormSup$'Worms')
  
  # CFU/5?l/Worm
  combiWormSup$'MeanW/Worm/25?l' <- (combiWormSup$'MeanW/Worm'* 10^combiWormSup$'Dilution.x')
  combiWormSup$'MeanS/Worm/25?l' <- (combiWormSup$'MeanS/Worm'* 10^combiWormSup$'Dilution.y')
  
  # CFU/100?l/Worm
  combiWormSup$'MeanW/Worm/100?l' <- (combiWormSup$'MeanW/Worm/25?l'*4)
  combiWormSup$'MeanS/Worm/100?l' <- (combiWormSup$'MeanS/Worm/25?l'*4)
  
  # CFUs
  combiWormSup$'CFUs' <- (combiWormSup$'MeanW/Worm/100?l'-combiWormSup$'MeanS/Worm/100?l')
  
  #final Tabel
  CFUs <- combiWormSup[,c(1:5,25)]
  
  Mut <- CFUs[!CFUs$MT == "dT", ]
  dT<- CFUs[CFUs$MT == "dT", ]
  
  # Merge worm and sup data frames
  totalCFU <- merge(Mut, dT, by = c("Replicate","Set","Worms","Sample"),
                    all.x = T)
  totalCFU<- totalCFU[!is.na(totalCFU$CFUs.y),]
  totalCFU$'CFUs'<-(totalCFU$'CFUs.x'+totalCFU$'CFUs.y')
  
  return(totalCFU)
}

#### apply function####

totalCFU <- CFUSCompColT(Competition)
totalCFU2 <- CFUSCompColT(Competition2)

totalCFU$Round<-as.factor(1)
totalCFU2$Round<-as.factor(2)
STPtotal<-rbind(totalCFU,totalCFU2)

#Fig. 2b
PlotCFUsT<- function(Data){
  Pop_Plot <- ggplot(Data,aes(x=Set,y=CFUs))+
    geom_boxplot(data=Data, aes(),lwd=0.5,fatten=1, colour= "black",outlier.shape = NA)+
    labs(fill="Morphotype")+
    #theme(legend.title=element_text(size=35),legend.text=element_text(size=35))+
    geom_point(data=Data,aes(linetype = "Replicate"),shape=21,size=2,fill="grey90", colour="black", position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0,seed=NA))+
    scale_y_log10(limits=c(1, 100000), breaks = c( 10, 100, 1000, 10000, 100000), labels = scientific_10 )+
    facet_grid(cols = vars(Round), scales = "free",space = "free")+
    my_theme()+
    #theme(text = element_text(size=35),panel.grid = element_line(color = "grey50", linetype = "dashed"),axis.text.x=element_text(angle=0,hjust=0.5))+ 
    labs(y="Short-term persistence [CFU/worm]",x="Bacteria")
  return(Pop_Plot)
}

#Fig. 2b
STPtotal$Set <- factor(STPtotal$Set,levels = c('MT48+dT','MT14+dT','algD+dT','alg44+dT','azu+dT','fimA+dT','fli+dT','usher+dT','pel+dT','fimA/usher+dT','fli/usher+dT','pel/usher+dT'),ordered = TRUE)
PlotCFUsT(STPtotal)

####Stats####

#Fig. 2a
data<-Pup
data$Set<-droplevels(data$Set)
data$Replicate<-droplevels(data$Replicate)
data$log_rel <- log10(data$Rel)

#Fig. 2a_2
data2<-Pup2
data2$Set<-droplevels(data2$Set)
data2$Replicate<-droplevels(data2$Replicate)
data2$log_rel <- log10(data2$Rel)

#define comparisons of interest
comp48 <-cbind.data.frame(rep('MT48+dT', 7), c('MT14+dT','algD+dT','alg44+dT','azu+dT','fimA+dT','fli+dT','usher+dT')) 
comp14 <-cbind.data.frame(rep('MT14+dT', 7), c('MT48+dT','algD+dT','alg44+dT','azu+dT','fimA+dT','fli+dT','usher+dT')) 

comp2_48 <-cbind.data.frame(rep('MT48+dT', 5), c('MT14+dT','fimA/usher+dT','fli/usher+dT','pel+dT','pel/usher+dT')) 
comp2_14 <-cbind.data.frame(rep('MT14+dT', 5), c('MT48+dT','fimA/usher+dT','fli/usher+dT','pel+dT','pel/usher+dT'))

stat<-function(comparisons,data,test){

  # Create object with test results
  res           <- as.data.frame(matrix(nrow = nrow(comparisons), ncol = 17))
  colnames(res) <- c("Group 1", "Group 2", "Shapiro-Wilk test statistic Group 1", "Shapiro-Wilk p value Group 1", "Shapiro-wilk test statistic Group 2", "Shapiro-Wilk p value Group 2", "F test: F statistic", "F test: num df", "F test: denom df", "F test: pvalue", "test for comparisons of central tendencies", "H0", "test statistic", "df", "p value", "FDR corrected p value", "reject H0")
  
  for (j in 1:nrow(res)) {
    # extract readings for group 1 and group 2
    data1       <- data$log_rel[data$Set == comparisons[j, 1] ]
    data2       <- data$log_rel[data$Set == comparisons[j, 2] ]
    
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
stat(comp48,data,"2a_anc")
stat(comp14,data,"2a_wspE")

S2a_anc<-read.csv("Stats_out/stats_2a_anc.csv",TRUE,",")
S2a_wspE<-read.csv("Stats_out/stats_2a_wspE.csv",TRUE,",")

stat(comp2_48,data2,"2a_anc_2")
stat(comp2_14,data2,"2a_wspE_2")

S2a_anc2<-read.csv("Stats_out/stats_2a_anc_2.csv",TRUE,",")
S2a_wspE2<-read.csv("Stats_out/stats_2a_wspE_2.csv",TRUE,",")

#Fig. S4, 2b
stat2<-function(comparisons,data,test){
  
  # Create object with test results
  res           <- as.data.frame(matrix(nrow = nrow(comparisons), ncol = 17))
  colnames(res) <- c("Group 1", "Group 2", "Shapiro-Wilk test statistic Group 1", "Shapiro-Wilk p value Group 1", "Shapiro-wilk test statistic Group 2", "Shapiro-Wilk p value Group 2", "F test: F statistic", "F test: num df", "F test: denom df", "F test: pvalue", "test for comparisons of central tendencies", "H0", "test statistic", "df", "p value", "FDR corrected p value", "reject H0")
  
  for (j in 1:nrow(res)) {
    # extract readings for group 1 and group 2
    data1       <- data$LOG_CFU[data$Set == comparisons[j, 1] ]
    data2       <- data$LOG_CFU[data$Set == comparisons[j, 2] ]
    
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

#Fig. S4
#define comparisons of interest
comp_solo<-cbind.data.frame(rep("MT48", 8), c('dT','MT14','algD','alg44','azu','fimA','fli','usher'))
comp2_solo<-cbind.data.frame(rep("MT48", 6), c('dT','MT14','pel','pel/usher','fimA/usher','fli/usher'))
comp_solo14<-cbind.data.frame(rep("MT14", 8), c('dT','MT48','algD','alg44','azu','fimA','fli','usher'))
comp2_solo14<-cbind.data.frame(rep("MT14", 6), c('MT48','dT','pel','pel/usher','fimA/usher','fli/usher'))

#execute function
stat2(comp_solo,CFUs_Comp,"S4_anc")
S4_anc<-read.csv("Stats_out/stats_S4_anc.csv",TRUE,",")

stat2(comp2_solo,CFUs_Comp2,"S4_anc_2")
S4_anc_2<-read.csv("Stats_out/stats_S4_anc_2.csv",TRUE,",")

stat2(comp_solo14,CFUs_Comp,"S4_wspE")
S4_wspE<-read.csv("Stats_out/stats_S4_wspE.csv",TRUE,",")

stat2(comp2_solo14,CFUs_Comp2,"S4_wspE_2")
S4_wspE2<-read.csv("Stats_out/stats_S4_wspE_2.csv",TRUE,",")

#Fig. 2b
# Pup
data<-totalCFU
data$Set<-droplevels(data$Set)
data$Replicate<-droplevels(data$Replicate)
data$LOG_CFU <- log10(data$CFUs)

# Pup2
data2<-totalCFU2
data2$Set<-droplevels(data2$Set)
data2$Replicate<-droplevels(data2$Replicate)
data2$LOG_CFU <- log10(data2$CFUs)

#execute function
stat2(comp48,data,"2b_anc")
stat2(comp14,data,"2b_wspE")

S2b_anc<-read.csv("Stats_out/stats_2b_anc.csv",TRUE,",")
S2b_wspE<-read.csv("Stats_out/stats_2b_wspE.csv",TRUE,",")

stat2(comp2_48,data2,"2b_anc_2")
stat2(comp2_14,data2,"2b_wspE_2")

S2b_anc2<-read.csv("Stats_out/stats_2b_anc_2.csv",TRUE,",")
S2b_wspE2<-read.csv("Stats_out/stats_2b_wspE_2.csv",TRUE,",")