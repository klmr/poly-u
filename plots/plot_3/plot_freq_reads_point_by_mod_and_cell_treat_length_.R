#Set working directory
setwd("/Users/marcos/Desktop/tailseq_celegans/plots/plot_3/data/")

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
                   aes(x = cell, y = frequency, fill = treatment)) +
    #geom_bar(stat = "identity", position=position_dodge(0.9), colour = "black", width = 0.8, size = thesize) +
    #geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9), size = thesize) +
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
  
  for(i in cells){
    topValue_1 <- max(data_short_p[data_short_p$cell == i,"mean_sd"])
    graph1 <- graph1 + geom_segment(x=j+0.78, xend=j+1.22, y=topValue_1 + constant, yend=topValue_1 + constant, size = thesize) + 
      geom_segment(x=j+0.78, xend=j+0.78, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) + 
      geom_segment(x=j+1.22, xend=j+1.22, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) +
      annotate("text",family="Arial", size=2.7, x=j+1,y=topValue_1 + constant,label= data_short_p[data_short_p$cell == i,"label"][1]) 
    j =j + 1 
  }
  return(graph1)
}

reads_all       <- read.csv("freq_reads_points_by_mod_and_cell_treat_p_change.csv")
reads_all$short <- rep("all", length(rownames(reads_all)))
reads_short     <- read.csv("freq_reads_points_by_mod_and_cell_treat_length_p_change.csv")
reads           <- full_join(reads_all, reads_short)
reads$short     <- as.factor(reads$short)  
reads <- reads[reads$mod %in% c("T", "TT", "TTT", "C", "G"),]

graph1 <- histogram_p(reads,6,1,"yes","T", c("liver", "bm", "mefs", "escs"))
graph1
graph2 <- histogram_p(reads,2.5,0.5,"no","T", c("liver", "bm", "mefs", "escs"))
graph2
graph3 <- histogram_p(reads,3,0.5,"all","T", c("liver", "bm", "mefs", "escs"))
graph3

graphs <- list(graph3, graph1, graph2)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_T.pdf", width=7, height=2.0, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,6,1,"yes","T", c("pachy", "liver", "bm", "mefs", "escs"))
graph1

graphs <- list(graph1)
graphs <- append(graphs, list(ncol=1))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_T_short_pachy.pdf", width=3.5, height=2.2, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,2,0.5,"yes","G", c("pachy", "liver", "bm", "mefs", "escs"))
graph1
graph2 <- histogram_p(reads,2,0.5,"no","G", c("pachy", "liver", "bm", "mefs", "escs"))
graph2
graph3 <- histogram_p(reads,2,0.5,"all","G", c("pachy", "liver", "bm", "mefs", "escs"))
graph3

graphs <- list(graph3, graph1, graph2)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_G.pdf", width=7.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,6,1,"yes","C", c("pachy", "liver", "bm", "mefs", "escs"))
graph1
graph2 <- histogram_p(reads,6,1,"no","C", c("pachy", "liver", "bm", "mefs", "escs"))
graph2
graph3 <- histogram_p(reads,6,1,"all","C", c("pachy", "liver", "bm", "mefs", "escs"))
graph3

graphs <- list(graph3, graph1, graph2)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_C.pdf", width=7.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,0.75,0.25,"yes","C", c("gv"))
graph1
graph2 <- histogram_p(reads,0.75,0.25,"no","C", c("gv"))
graph2
graph3 <- histogram_p(reads,0.75,0.25,"all","C", c("gv"))
graph3

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_C_gv.pdf", width=2.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()


graph1 <- histogram_p(reads,1,0.25,"yes","G", c("gv"))
graph1
graph2 <- histogram_p(reads,1,0.25,"no","G", c("gv"))
graph2
graph3 <- histogram_p(reads,1,0.25,"all","G", c("gv"))
graph3

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_G_gv.pdf", width=2.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,2,0.5,"yes","T", c("gv"))
graph1
graph2 <- histogram_p(reads,2,0.5,"no","T", c("gv"))
graph2
graph3 <- histogram_p(reads,2,0.5,"all","T", c("gv"))
graph3

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_T_gv.pdf", width=2.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- histogram_p(reads,1.25,0.25,"yes","TT", c("gv"))
graph1
graph2 <- histogram_p(reads,1.25,0.25,"no","TT", c("gv"))
graph2
graph3 <- histogram_p(reads,1.25,0.25,"all","TT", c("gv"))
graph3

graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)
pdf("../results/freq_reads_point_by_mod_and_cell_treat_length_TT_gv.pdf", width=2.5, height=2.5, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

