## -----------------------------------------------------------------------------
## 
## Analysis of the International Trade Network data
##
## This script plots the clustering coefficient for the original data
## for each year together with 1 standard deviation confidence
## intervals for the clustering coefficient from surrogates
## constructed using (i) random shuffling of edge weights and where
## (ii) vertex strenghts are preserved exactly and (iii) vertex
## strengths are allowed to vary on an interval.
##
## The data needed for the plots is also saved as an rds-file.
## 
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla analyse_itn_results.R <data_dir> <results_dir> <figure_dir>
##
## data_dir    : path to the directory holding the ITN trade data
## results_dir : path to the directory holding the results of the ITN analyses
## figure_dir  : path to the directory where the result figures should be saved
## -----------------------------------------------------------------------------

source("utilities.R")
library(latex2exp)
library(ggplot2)

## -----------------------------------------------------------------------------
## Get command-line arguments
## -----------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
    stop("Not enough arguments")

data_dir    <- args[1]
results_dir <- args[2]
figure_dir  <- args[3]

## -----------------------------------------------------------------------------
## Set the random seed
## -----------------------------------------------------------------------------
set.seed(42)

## -----------------------------------------------------------------------------
## Calculate the clustering coefficient in the original data
## and for weight-preserving surrogates (shuffling edge weights)
## -----------------------------------------------------------------------------

data_itn <- readRDS(file.path(data_dir, "itn_data.rds"))

cc_orig_undir <- sapply(names(data_itn$itn_u), function(y) get_cc(data_itn$itn_u[[y]], max_norm = TRUE, directed = FALSE))
cc_orig_dir   <- sapply(names(data_itn$itn_d), function(y) get_cc(data_itn$itn_d[[y]], max_norm = TRUE, directed = TRUE))

## (1) get the nimer of surrogates needed
N_s <- length(readRDS(file = list.files(results_dir, full.names = TRUE)[1])$cc)

## shuffle edge weights
cc_shuffle_undir <- t(sapply(names(data_itn$itn_u), function(y) calc_stats(replicate(N_s, get_cc(shuffle_weights(data_itn$itn_u[[y]]), max_norm = TRUE, directed = FALSE)), y)))
cc_shuffle_dir   <- t(sapply(names(data_itn$itn_d), function(y) calc_stats(replicate(N_s, get_cc(shuffle_weights(data_itn$itn_d[[y]]), max_norm = TRUE, directed = TRUE)), y)))

## -----------------------------------------------------------------------------
## Read the results
## -----------------------------------------------------------------------------

## There are four types of data: {dir, undir} x {fixed, interval}. Read each type.
cc_undir_fixed    <- read_itn_results(basepath = results_dir, network_type = "undir", sampler_type = "fixed")
cc_undir_interval <- read_itn_results(basepath = results_dir, network_type = "undir", sampler_type = "interval")

cc_dir_fixed      <- read_itn_results(basepath = results_dir, network_type = "dir", sampler_type = "fixed")
cc_dir_interval   <- read_itn_results(basepath = results_dir, network_type = "dir", sampler_type = "interval")

## -----------------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------------
saveRDS(list("cc_orig_undir"     = cc_orig_undir,
             "cc_orig_dir"       = cc_orig_dir,
             "cc_shuffle_undir"  = cc_shuffle_undir,
             "cc_shuffle_dir"    = cc_shuffle_dir,
             "cc_undir_fixed"    = cc_undir_fixed,
             "cc_undir_interval" = cc_undir_interval,
             "cc_dir_fixed"      = cc_dir_fixed,
             "cc_dir_interval"   = cc_dir_interval),
        file = file.path(figure_dir, "itn_results.rds"), compress = "xz")

## -----------------------------------------------------------------------------
## Plot results
## -----------------------------------------------------------------------------

ylim_undir <- range(c(cc_shuffle_undir[, 2], cc_undir_fixed[, 3], cc_undir_interval[, 3]))
ylim_dir   <- range(c(cc_shuffle_dir[, 2], cc_dir_fixed[, 3], cc_dir_interval[, 3]))

xlab <- "year"
ylab <- TeX("$\\bar{C}$")

p1 <- make_itn_plot(cc_orig_undir, cc_shuffle_undir, cc_undir_fixed, cc_undir_interval, ylim = ylim_undir, main = "undir - fixed", xlab = xlab, ylab = ylab, notitle = TRUE)
p2 <- make_itn_plot(cc_orig_dir, cc_shuffle_dir, cc_dir_fixed, cc_dir_interval, ylim = ylim_dir, main = "dir - fixed", xlab = xlab, ylab = ylab, notitle = TRUE)

ggsave(plot = p1, filename = file.path(figure_dir, "cc_undir.pdf"), width = 150, height = 100, units = "mm")
ggsave(plot = p2, filename = file.path(figure_dir, "cc_dir.pdf"), width = 150, height = 100, units = "mm")
