#Set working directory
setwd("/Users/Marcos/Dropbox/Lab/shared/collaborations/tailseq_celegans/poly-u/plots/")
reads <- read.csv("../../plots/sample_table.csv")

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

#reads  <- readsx
readsx <- reads
reads$palen  <- ifelse(readsx$cell == "liver" | readsx$cell == "escs", readsx$palen + 1, readsx$palen) 

doGraph <- function(cells, colors){
  grays <- colors
  thesize = 0.468
  reads  <- reads %>% mutate(cell = factor(cell, levels = cells))
  graph1 <- ggplot(aes(x = palen, y = 100 * ..density..), 
                  data = subset(reads, palen < 100 & treatment == "ctrl" & cell %in% cells)) +
    geom_freqpoly(aes(color = cell), binwidth=6) + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 80, 20), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0, 1.5, 0.3), expand = c(0, 0)) +
    geom_vline(xintercept = 30, linetype = "longdash") +
    scale_color_manual(values=grays) +
    theme_classic(base_size = 10) + #, text=element_text(family="Arial", size=8)) +
    coord_cartesian(ylim = c(0, 1.5), xlim = c(0,80)) +
    xlab('Poly(A) length (nt)') + ylab("Density of Reads [AU]") +
    theme(legend.position = "right", legend.title = element_blank(),
          text=element_text(family="Arial", size=8), 
          axis.line.x = element_line(size = thesize), axis.line.y = element_line(size = thesize)) 
  return(graph1)
}

graph1 <- doGraph(c("liver", "bm", "mefs", "escs"), gray(c(0.2,0.4,0.6,0.8)))
graph1 
ggsave("../../plots/plot_1/results/density_read_by_length_and_cell_not_pachy.pdf", graph1, width=2.5, height=2.0)
rm(graph1)

graph1 <- doGraph(c("pachy","liver", "bm", "mefs", "escs"), c("red", gray(c(0.2,0.4,0.6,0.8))))
graph1 
ggsave("../../plots/plot_1/results/density_read_by_length_and_cell_pachy.pdf", graph1, width=2.8, height=2.2)
rm(graph1)

graph1 <- doGraph(c("gv"), c("red"))
graph1 
ggsave("../../plots/plot_1/results/density_read_by_length_and_cell_gv.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraph(c("gv", "pachy","liver", "bm", "mefs", "escs"), c("red", "green4", gray(c(0.2,0.3,0.4,0.6,0.8))))
graph1 
ggsave("../../plots/plot_1/results/density_read_by_length_and_cell_all.pdf", graph1, width=3.2, height=2.7)
rm(graph1)


