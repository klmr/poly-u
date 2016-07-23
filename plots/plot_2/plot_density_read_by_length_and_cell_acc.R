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

#plot_density_read_by_length_and_cell_acc


readsx <- reads
reads$palen  <- ifelse(readsx$cell == "liver" | readsx$cell == "escs", readsx$palen + 1, readsx$palen) 
reads  <- readsx

doGraph <- function(a_cell, a_select, a_treat, scale, breack){

  thesize = 0.234
  reads$select <- reads[,c(a_select)]

  graph1 <- ggplot(aes(x = palen, y = 100 * ..density..), 
                   data = subset(reads, palen < 100 & treatment == a_treat & cell == a_cell)) +
    geom_freqpoly(aes(color = select), binwidth=6) + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 80, 20), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0, scale, breack), expand = c(0, 0)) +
    geom_vline(xintercept = 30, linetype = "longdash") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic(base_size = 10) +
    coord_cartesian(ylim = c(0, scale), xlim = c(0,80)) +
    theme(legend.position = "bottom", 
          text=element_text(family="Arial", size=8), 
          axis.line.x = element_line(size = thesize), axis.line.y = element_line(size = thesize)) +
    xlab('Poly(A) length (nt)') + ylab("Density of Reads [AU]")

  return(graph1)
}

graph1 <- doGraph("pachy", "pachy", "ctrl", 1.5, 0.3)
graph1 
ggsave("../../plots/plot_2/results/plot_poly_read_by_length_and_cell_acc_pachy_ctrl.pdf", graph1, width=2.0, height=2.8)
rm(graph1)

graph1 <- doGraph("pachy", "pachy", "dcko", 1.5, 0.3)
graph1 
ggsave("../../plots/plot_2/results/plot_poly_read_by_length_and_cell_acc_pachy_dcko.pdf", graph1, width=2.0, height=2.8)
rm(graph1)

graph1 <- doGraph("gv", "gv", "ctrl", 2.75, 0.25)
graph1 
ggsave("../../plots/plot_2/results/plot_poly_read_by_length_and_cell_acc_gv_ctrl.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraph("gv", "gv", "dcko", 2.75, 0.25)
graph1 
ggsave("../../plots/plot_2/results/plot_poly_read_by_length_and_cell_acc_gv_dcko.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

#####

doGraphDensity <- function(a_cell, a_select, a_treat, scale, breack){
  
  thesize = 0.234
  reads$select <- reads[,c(a_select)]
  
  graph1 <- ggplot(aes(x = palen, y = 100 * ..density..), 
                   data = subset(reads, palen < 80 & treatment %in% a_treat & cell == a_cell)) +
    geom_density(aes(color = select)) + # , binwidth=6) + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 80, 20), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0, scale, breack), expand = c(0, 0)) +
    geom_vline(xintercept = 30, linetype = "longdash") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic(base_size = 10) +
    coord_cartesian(ylim = c(0, scale), xlim = c(0,80)) +
    theme(legend.position = "bottom", 
          text=element_text(family="Arial", size=8), 
          axis.line.x = element_line(size = thesize), axis.line.y = element_line(size = thesize)) +
    xlab('Poly(A) length (nt)') + ylab("Density of Reads [AU]")
  
  return(graph1)
}

data = subset(reads, palen < 80 & treatment %in% c("dcko") & cell == "gv")

graph1 <- doGraphDensity("gv", "gv", c("ctrl"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_read_by_length_and_cell_acc_gv_ctrl.pdf", graph1, width=2.3, height=2.7)
rm(graph1)

graph1 <- doGraphDensity("gv", "gv", c("dcko"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_read_by_length_and_cell_acc_gv_dcko.pdf", graph1, width=2.3, height=2.7)
rm(graph1)

graph1 <- doGraphDensity("pachy", "pachy", c("ctrl"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_read_by_length_and_cell_acc_pachy_ctrl.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensity("pachy", "pachy", c("dcko"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_read_by_length_and_cell_acc_pachy_dcko.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensity("pachy", "treatment", c("dcko","ctrl"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_read_by_length_and_cell_treatment_pachy.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

#####

doGraphDensityAll <- function(a_cell, a_select, a_treat, scale, breack){
  
  thesize = 0.234
  reads$select <- reads[,c(a_select)]
  
  graph1 <- ggplot(aes(x = palen, y = 100 * ..density..), 
                   data = subset(reads, palen < 100 & treatment %in% a_treat & cell == a_cell)) +
    geom_density(aes(color = select)) + # , binwidth=6) + 
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0, scale, breack), expand = c(0, 0)) +
    geom_vline(xintercept = 30, linetype = "longdash") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic(base_size = 10) +
    coord_cartesian(ylim = c(0, scale), xlim = c(0,100)) +
    theme(legend.position = "bottom", 
          text=element_text(family="Arial", size=8), 
          axis.line.x = element_line(size = thesize), axis.line.y = element_line(size = thesize)) +
    xlab('Poly(A) length (nt)') + ylab("Density of Reads [AU]")
  
  return(graph1)
}

graph1 <- doGraphDensityAll("gv", "gv", c("ctrl"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_all_read_by_length_and_cell_acc_gv_ctrl.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensityAll("gv", "gv", c("dcko"), 3.5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_all_read_by_length_and_cell_acc_gv_dcko.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensityAll("pachy", "pachy", c("ctrl"), 10, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_all_read_by_length_and_cell_acc_pachy_ctrl.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensityAll("pachy", "pachy", c("dcko"), 5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_all_read_by_length_and_cell_acc_pachy_dcko.pdf", graph1, width=2.7, height=2.7)
rm(graph1)

graph1 <- doGraphDensityAll("pachy", "treatment", c("dcko","ctrl"), 5, 0.5)
graph1 
ggsave("../../plots/plot_2/results/plot_density_all_read_by_length_and_cell_treatment_pachy.pdf", graph1, width=2.7, height=2.7)
rm(graph1)
