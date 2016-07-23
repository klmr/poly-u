#Set working directory

setwd("/Users/marcos/Desktop/tailseq_celegans/plots/plot_5/data/")

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

reads <- read.csv("freq_points_by_mod_and_cell_treat_group_p_change.csv")
reads <- reads[, !colnames(reads) %in% c("X", "X.1") ]
reads <- mutate(reads, type = paste(Treatment, Tissue, Group, sep = "_")) 
reads <- mutate(reads, type = factor(type, levels = type))
reads <- reads %>% mutate(label = paste(as.character(signif(change, digits = 2)), "X" , "\n", star, sep =""))
reads <- reads %>% mutate(mean_sd = mean + sd)



do_plot <- function(cells, mods, topValue, breaks, the_colors){

  reads <- mutate(reads, Modification = factor(Modification, levels = mods))
  
  constant <- 1
  thesize = 0.234 

  graph1 <- ggplot(subset(reads, Modification %in% mods & Tissue %in% cells), 
                   aes(x = factor(type), y = Frequency, fill = factor(Modification))) +
    geom_dotplot(binaxis = "y", stackdir = "center") +
    stat_summary_bin(aes(y = Frequency, colors = "black"), fun.y = "mean", geom = "bar", width = 0.2) +
    stat_summary(fun.data=mean_sdl,  fun.args = list(mult = 1), color = "black", size = 0.5) +
  
    scale_y_continuous(breaks=seq(0,topValue,breaks), expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    coord_cartesian(ylim = c(0, topValue)) +
    scale_fill_manual(values=the_colors, guide = FALSE) +
    theme_classic() + 
    ggtitle(paste(cells[1], mods, sep=" ")) + 
    ylab("Frequency of Reads [%]") + 
    theme(strip.background = element_blank(), text=element_text(family="Arial", size=8), plot.margin=unit(c(0.15,0.15,0.15,0.15), "cm"),
          legend.position = c(0.9,0.9),  legend.title=element_blank(),
          axis.line = element_line(size = thesize),
          axis.ticks.x=element_blank(),  axis.title.x = element_blank(), 
          axis.ticks.y=element_line(size = thesize),
          axis.line.y = element_line(size = thesize),
          axis.line.x = element_line(size = thesize))  

    for(k in mods){
      j <- 0
      reads_short <- subset(reads, Modification == k & Tissue %in% cells)
      for(i in c("All", "let-7")){
        topValue_1 <- max(reads_short[reads_short$Group == i,"mean_sd"])
        graph1 <- graph1 + geom_segment(x=j+1, xend=j+2, y=topValue_1 + constant, yend=topValue_1 + constant, size = thesize) +
          geom_segment(x=j+1, xend=j+1, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) +
          geom_segment(x=j+2, xend=j+2, y=topValue_1 + constant, yend=topValue_1 + constant/2, size = thesize) +
          annotate("text",family="Arial", size=2.7, x=j+1.5,y=topValue_1 + constant,label= reads_short[reads_short$Group == i,"label"][1])
        j =j + 2
      }
    }
  
  return(graph1)
}

mods   <- c("U_all", "PolyU")
the_colors <- c("slategray1", "steelblue3")
graph1 <- do_plot(c("Liver"), mods, 12, 1, the_colors)  
graph2 <- do_plot(c("BM"), mods, 10, 1, the_colors)  
graph3 <- do_plot(c("MEFs"), mods, 13, 1, the_colors)  
graph4 <- do_plot(c("ESCs"), mods, 10, 1, the_colors)  
graphs <- list()
graphs <- list(graph1, graph2, graph3, graph4)
graphs <- append(graphs, list(ncol=4))
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_u.pdf", width=7, height=1.85, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- do_plot(c("Liver"), c("A_all", "PolyA"), 13, 1, c("peachpuff","sandybrown"))  
graph2 <- do_plot(c("Liver"), c("C_all", "PolyC"), 26, 2, c("darkseagreen2","palegreen4"))  
graph3 <- do_plot(c("Liver"), c("G_all", "PolyG"), 7, 1, c("thistle","thistle4"))  
graphs <- list()
graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_Liver_not_u.pdf", width=5.2, height=1.85, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- do_plot(c("BM"), c("A_all", "PolyA"), 12, 1, c("peachpuff","sandybrown"))  
graph2 <- do_plot(c("BM"), c("C_all", "PolyC"), 16, 2, c("darkseagreen2","palegreen4"))  
graph3 <- do_plot(c("BM"), c("G_all", "PolyG"), 5, 1, c("thistle","thistle4"))  
graphs <- list()
graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_BM_not_u.pdf", width=5.2, height=1.85, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- do_plot(c("MEFs"), c("A_all", "PolyA"), 18, 2, c("peachpuff","sandybrown"))  
graph2 <- do_plot(c("MEFs"), c("C_all", "PolyC"), 20, 2, c("darkseagreen2","palegreen4"))  
graph3 <- do_plot(c("MEFs"), c("G_all", "PolyG"), 6, 1, c("thistle","thistle4"))  
graphs <- list()
graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_MEFs_not_u.pdf", width=5.2, height=1.85, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()

graph1 <- do_plot(c("ESCs"), c("A_all", "PolyA"), 12, 1, c("peachpuff","sandybrown"))  
graph2 <- do_plot(c("ESCs"), c("C_all", "PolyC"), 28, 4, c("darkseagreen2","palegreen4"))  
graph3 <- do_plot(c("ESCs"), c("G_all", "PolyG"), 16, 4, c("thistle","thistle4"))  
graphs <- list()
graphs <- list(graph1, graph2, graph3)
graphs <- append(graphs, list(ncol=3))
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_ESCs_not_u.pdf", width=5.2, height=1.85, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()


graph1 <- do_plot(c("Pachy"), c("A_all", "PolyA"), 12, 1, c("peachpuff","sandybrown"))  
graph2 <- do_plot(c("Pachy"), c("C_all", "PolyC"), 30, 5, c("darkseagreen2","palegreen4"))  
graph3 <- do_plot(c("Pachy"), c("G_all", "PolyG"), 12, 3, c("thistle","thistle4"))  
graph4 <- do_plot(c("Pachy"), c("U_all", "PolyU"), 25, 5, c("slategray1", "steelblue3")) 
graphs <- list()
graphs <- list(graph1, graph2, graph3, graph4)
graphs <- append(graphs, list(ncol=2, nrow=2)
do.call(grid.arrange, graphs)

do.call(grid.arrange, graphs)
pdf("../results/freq_mod_points_by_cell_and_type_treat_Pachy_all.pdf", width=3.4, height=4.2, useDingbats=FALSE)
do.call(grid.arrange, graphs)
dev.off()





