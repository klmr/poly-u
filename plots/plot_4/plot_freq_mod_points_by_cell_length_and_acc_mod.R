#Set working directory
setwd("/Volumes/Promise_Pegasus/marcos_dropbox/Dropbox/Lab/Marcos_files/Scripts_Data/scripts/tail_seq/cambridge_first_week_2016/TAIL-Seq_graphs/all/data/")

setwd("/Users/marcos/Desktop/tailseq_celegans/plots/plot_4/data/")

library(dplyr)                                   #load libraries for data wrangling and ploting
library(tidyr)
#library(plyr)
library(ggplot2)
library(gridExtra)
library(extrafont)                               #Loads fonts to include Arial
library(extrafontdb)
font_import(recursive = FALSE)                   #NOTE: Needs to stop at this point and type y in the console
loadfonts(device = "postscript")
loadfonts()

#Main function to plot
histogram_p <- function(data_short_p, topValue, breaks, tails, modif, cells){
  
  data_short_p <- reads
  constant <- topValue/10
  j=0
  thesize = 0.234 #
  
  data_short_p <- data_short_p %>% mutate(label = paste(as.character(signif(change, digits = 2)), "X" , "\n", star, sep =""))
  data_short_p <- data_short_p %>% mutate(mean_sd = mean + sd)
  data_short_p <- subset(data_short_p, mod %in% c(modif) & short %in% tails & cell %in% cells) 
  data_short_p <- data_short_p %>% mutate(cell = factor(cell, levels = cells))
  
  graph1 <- ggplot(data_short_p, 
                   aes(x = pachy, y = frequency, fill = treatment)) +
    geom_dotplot(binaxis = "y", stackdir = "center", position=position_dodge(width=0.9)) +
    stat_summary_bin(aes(y = frequency, colors = "black"), fun.y = "mean", geom = "bar", width = 0.2, position = position_dodge(width=0.9)) +
    stat_summary(fun.data=mean_sdl,  fun.args = list(mult = 1), position=position_dodge(width=0.9), color = "black", size = 0.5) +
    
    scale_y_continuous(breaks=seq(0,topValue,breaks), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    coord_cartesian(ylim = c(0, topValue)) +
    scale_fill_manual(values=c("white", "black")) +
    theme_classic() + 
    ggtitle(paste(tails, modif, sep=" ")) + 
    ylab("Frequency of Reads [%]") + 
    theme(strip.background = element_blank(), text=element_text(family="Arial", size=8), plot.margin=unit(c(0.15,0.15,0.15,0.15), "cm"),
          legend.position = c(0.9,0.9),  legend.title=element_blank(),
          axis.line = element_line(size = thesize),
          axis.ticks.x=element_blank(),  axis.title.x = element_blank(), 
          axis.ticks.y=element_line(size = thesize),
          axis.line.y = element_line(size = thesize),
          axis.line.x = element_line(size = thesize))  
  
  for(i in c("no", "yes")){
    topValue_1 <- max(data_short_p[data_short_p$pachy == i,"mean_sd"])
    graph1 <- graph1 + geom_segment(x=j+0.78, xend=j+1.22, y=topValue_1 + constant, yend=topValue_1 + constant, size = thesize) +
      geom_segment(x=j+0.78, xend=j+0.78, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) +
      geom_segment(x=j+1.22, xend=j+1.22, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) +
      annotate("text",family="Arial", size=2.7, x=j+1,y=topValue_1 + constant,label= data_short_p[data_short_p$pachy == i,"label"][1])
    j =j + 1
  }
  return(graph1)
}

reads_all       <- read.csv("freq_reads_points_by_mod_and_cell_treat_acc_p_change_x_.csv")
reads_all$short <- rep("all", length(rownames(reads_all)))
reads_short     <- read.csv("freq_reads_points_by_mod_and_cell_treat_acc_length_p_change_x_.csv")
reads           <- full_join(reads_all, reads_short)
reads$short     <- as.factor(reads$short)  
reads <- reads[reads$mod %in% c("T", "TT", "TTT", "C", "G"),]

graph1 <- histogram_p(reads,8,1,"yes","T", c("pachy"))
graph2 <- histogram_p(reads,8,1,"no","T", c("pachy"))

graphs <- list(graph1, graph2)
graphs <- append(graphs, list(ncol=2))
do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_length_and_acc_mod_T_not_all_pachy.pdf", width=3.5, height=2.2, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,8,1,"yes","T", c("pachy"))
graph2 <- histogram_p(reads,8,1,"no","T", c("pachy"))
graph3 <- histogram_p(reads,8,1,"all","T", c("pachy"))

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_length_and_acc_mod_T_pachy.pdf", width=7, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()


graph1 <- histogram_p(reads,8,1,"yes","T", c("liver"))
graph2 <- histogram_p(reads,8,1,"no","T", c("liver"))
graph3 <- histogram_p(reads,8,1,"all","T", c("liver"))

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_length_and_acc_mod_T_liver.pdf", width=7, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()


