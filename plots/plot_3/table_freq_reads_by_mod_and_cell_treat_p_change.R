library(dplyr)
library(tidyr)

setwd("/Users/Marcos/Dropbox/Lab/shared/collaborations/tailseq_celegans/poly-u/plots/")
reads <- read.csv("../../plots/sample_table.csv")


reads           <- reads_filt_mrna[, c("mod","cell","treatment","sample")]
rm(reads_filt_mrna_sample)

reads_all       <- reads[,c("cell","treatment","sample")]
reads_all       <- summarise(group_by(reads_all, cell, treatment, sample), all = n())
reads           <- summarise(group_by(reads, cell, treatment, sample, mod), mod_count = n())
reads_summary   <- full_join(reads_all, reads)
rm(reads_all, reads)

reads_summary   <- mutate(reads_summary, frequency = 100*mod_count/all)
reads_summary   <- reads_summary[,c("cell", "treatment", "sample", "mod", "frequency")]
reads_stats     <- summarise(group_by(reads_summary, cell, treatment, mod), mean = mean(frequency), sd = sd(frequency), frequencies=list(frequency))

#fold change
{
reads_stats_f   <- reads_stats
reads_stats_f   <- reads_stats_f[,c("cell", "treatment", "mod", "mean")]
reads_stats_f   <- spread(reads_stats_f, treatment, mean)
  
vector <- c()
for( i in 1:length(reads_stats_f$ctrl)){
  x <- as.vector(unlist(reads_stats_f[i,"ctrl"]))
  y <- as.vector(unlist(reads_stats_f[i,"dcko"]))
  change <- ifelse(!is.na(x)  & !is.na(y), 
                    ifelse((x/y) > 1, x/y, -y/x)
                    , NA)
  vector <- c(vector, change)
}
reads_stats_f$change <- vector
rm(x, y, vector, change)
reads_stats_f <- reads_stats_f[, c("cell", "mod", "change")]
reads_stats   <- full_join(reads_stats, reads_stats_f)
}

#t-test
{
reads_stats_t   <- reads_stats
reads_stats_t   <- reads_stats_t[,c("cell", "treatment", "mod", "frequencies")]
reads_stats_t   <- spread(reads_stats_t, treatment, frequencies)

vector <- c()
for( i in 1:length(reads_stats_t$ctrl)){
  x <- as.vector(unlist(reads_stats_t[i,"ctrl"]))
  y <- as.vector(unlist(reads_stats_t[i,"dcko"]))
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
reads_stats_t <- reads_stats_t[, c("cell", "mod", "pval", "star")]

reads_stats   <- full_join(reads_stats, reads_stats_t)
reads_stats   <- reads_stats[,!(names(reads_stats) %in% c("frequencies"))]
}

reads_summary <- full_join(reads_summary, reads_stats, by= c("cell", "treatment", "mod"))

write.csv(reads_summary, file = "../../plots/plot_3/data/freq_reads_points_by_mod_and_cell_treat_p_change.csv")

