## -----------------------------------------------------------------------------
## Analyse the convergence results
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla analyse_results.R datapath1 datapath2 outputpath
##         
## -----------------------------------------------------------------------------

source("init.R")
source("utilities.R")
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(pracma)
library(xtable)

## --------------------------------------------------------------------------------
## Get command-line arguments
## --------------------------------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
respath    <- args[1]                            ## directory with convergence results
respath2   <- args[2]                            ## directory with running time results
outputpath <- args[3]                            ## directory where to store the results

## --------------------------------------------------------------------------------
## Plot the results separately for experiment 1 and experiment 2
## --------------------------------------------------------------------------------

p1 <- load_results_and_plot(basepath = respath, fpattern = "_convergence.rds", log_ds = TRUE)
p2 <- load_results_and_plot(basepath = respath, fpattern = "_convergence_nwfix.rds", log_ds = TRUE)

ggsave(file = file.path(outputpath, "convergence.pdf"), plot = p1, width = 130, height = 100, units = "mm")
ggsave(file = file.path(outputpath, "convergence_nwfix.pdf"), plot = p2, width = 130, height = 100, units = "mm")

## --------------------------------------------------------------------------------
## Make a table with dataset properties
## --------------------------------------------------------------------------------

res <- make_network_property_table(respath = respath2)
print(xtable(res, display = c("s", rep("e", 6), rep("f", 4))), booktabs = TRUE,  math.style.exponents = TRUE, sanitize.rownames.function = bold)

## --------------------------------------------------------------------------------
