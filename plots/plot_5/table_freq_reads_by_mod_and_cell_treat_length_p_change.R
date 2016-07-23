library(dplyr)
library(tidyr)

setwd("/Users/marcos/Desktop/tailseq_celegans/plots/plot_5/data/")

reads           <- read.csv("All miRNA modifications.csv", header = TRUE, row.names=NULL, sep = ";")
reads$Frequency <- as.numeric(gsub(",", ".", as.vector(reads$Frequency)))
reads           <- subset(reads, Treatment %in% c("Ctrl", "cKO"))
reads           <- reads[,!colnames(reads) %in% c("X")]

all <- ungroup(summarise(group_by(reads, Nucleotide, Tissue, Treatment, Group, Sample), Frequency = sum(Frequency)))
all <- mutate(all, Modification = paste(Nucleotide, "all", sep = "_"))

reads <- bind_rows(reads, all)

reads_stats     <- summarise(group_by(reads, Group, Tissue, Treatment, Modification), mean = mean(Frequency), sd = sd(Frequency), frequencies=list(Frequency))

#fold change
{
  reads_stats_f   <- reads_stats
  reads_stats_f   <- reads_stats_f[,c("Tissue", "Treatment", "Group", "Modification", "mean")]
  reads_stats_f   <- spread(reads_stats_f, Treatment, mean)
  
  vector <- c()
  for( i in 1:length(reads_stats_f$Ctrl)){
    x <- as.vector(unlist(reads_stats_f[i,"Ctrl"]))
    y <- as.vector(unlist(reads_stats_f[i,"cKO"]))
    change <- ifelse(!is.na(x)  & !is.na(y), 
                     ifelse((x/y) > 1, x/y, -y/x)
                     , NA)
    vector <- c(vector, change)
  }
  reads_stats_f$change <- vector
  rm(x, y, vector, change)
  reads_stats_f <- reads_stats_f[, c("Tissue", "Group", "Modification", "change")]
  reads_stats   <- full_join(reads_stats, reads_stats_f)
}

#t-test
{
  reads_stats_t   <- reads_stats
  reads_stats_t   <- reads_stats_t[,c("Tissue", "Treatment", "Group", "Modification", "frequencies")]
  reads_stats_t   <- spread(reads_stats_t, Treatment, frequencies)
  
  vector <- c()
  for( i in 1:length(reads_stats_t$Ctrl)){
    x <- as.vector(unlist(reads_stats_t[i,"Ctrl"]))
    y <- as.vector(unlist(reads_stats_t[i,"cKO"]))
    p_value <- ifelse(length(x) > 1  & length(y) > 1, t.test(x, y, var.equal=TRUE)$p.value, NA)
    vector <- c(vector, p_value)
  }
  reads_stats_t$pval <- vector
  rm(x, y, vector, p_value)
  
  sig <- function(i){
    res <- "NS"
    if(i < 0.05) res = "*"  
    if(i < 0.01) res = "**" 
    if(i < 0.001) res = "***" 
    return( res)
  }
  vector <- c()
  for(i in 1:length(reads_stats_t$pval)){
    star <- ifelse(!is.na(reads_stats_t$pval[i]), sig(reads_stats_t$pval[i]),NA)
    vector <- c(vector, star)
  }
  reads_stats_t$star <- vector
  rm(vector)
  reads_stats_t <- reads_stats_t[, c("Tissue", "Group", "Modification", "pval", "star")]
  
  reads_stats   <- full_join(reads_stats, reads_stats_t)
  reads_stats   <- reads_stats[,!(names(reads_stats) %in% c("frequencies"))]
}

reads_summary <- full_join(reads, reads_stats, by= c("Tissue", "Treatment", "Modification", "Group"))

write.csv(reads_summary, file = "freq_points_by_mod_and_cell_treat_group_p_change.csv")


