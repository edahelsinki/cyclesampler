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
library(latex2exp)
library(cowplot)

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

dslist <- c("Last.fm", "MovieLens_100k", "FineFoods", "MovieLens_1M", "BookCrossing", "MovieLens_20M", "TasteProfile")

plist_1 <- lapply(dslist, function(i) load_results_and_plot_l2_norm(basepath = respath, dsname = i, suffix = "convergence", log_ds = TRUE))
plist_2 <- lapply(dslist, function(i) load_results_and_plot_l2_norm(basepath = respath, dsname = i, suffix = "convergence_nwfix", log_ds = TRUE))

p1 <- do.call("plot_grid", c(plist_1, ncol = 2))
p2 <- do.call("plot_grid", c(plist_2, ncol = 2))

save_plot(filename = file.path(outputpath, "convergence.png"), plot = p1, ncol = 2, base_width = 4.2, base_height = 8, dpi = 300)
save_plot(filename = file.path(outputpath, "convergence_nwfix.png"), plot = p2, ncol = 2, base_width = 4.2, base_height = 8, dpi = 300)

## --------------------------------------------------------------------------------
## Make a table with dataset properties
## --------------------------------------------------------------------------------

res <- make_network_property_table(respath = respath2)
print(xtable(res, display = c("s", rep("e", 6), rep("f", 4))), booktabs = TRUE,  math.style.exponents = TRUE, sanitize.rownames.function = bold)

## --------------------------------------------------------------------------------
