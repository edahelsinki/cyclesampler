## -----------------------------------------------------------------------------
##
## Analysis of the International Trade Network data
##
## This script calculates the clustering coefficient for the original data
## for each year together with 1 standard deviation confidence
## intervals for the clustering coefficient from surrogates
## constructed using (i) random shuffling of edge weights and where
## (ii) vertex strenghts are preserved exactly and (iii) vertex
## strengths are allowed to vary on an interval.
##
## This script generates data needed for plots. The plots can then be generated using
## the script plot_itn_results.R.
##
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla analyse_itn_results.R <normalise_cc> <data_dir> <results_dir> <output_dir>
##
## normalise_cc : Boolean (0 or 1) indicating whether the clustering coefficient should be normalised or not
## data_dir     : path to the directory holding the ITN trade data
## results_dir  : path to the directory holding the results of the ITN analyses
## output_dir   : path to the directory where the results (rds-files) should be saved
## -----------------------------------------------------------------------------

source("utilities.R")

## -----------------------------------------------------------------------------
## Get command-line arguments
## -----------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
    stop("Not enough arguments")

normalise_cc <- as.logical(as.numeric(args[1])) ## normalise clustering coefficient using maximum edge weight (Boolean)
data_dir     <- args[2]
results_dir  <- args[3]
output_dir   <- args[4]

## -----------------------------------------------------------------------------
## Set the random seed
## -----------------------------------------------------------------------------
set.seed(42)

## -----------------------------------------------------------------------------
## Calculate the clustering coefficient in the original data
## and for weight-preserving surrogates (shuffling edge weights)
## -----------------------------------------------------------------------------

data_itn <- readRDS(file.path(data_dir, "itn_data.rds"))

cc_orig_undir <- sapply(names(data_itn$itn_u), function(y) get_cc(data_itn$itn_u[[y]], max_norm = normalise_cc, directed = FALSE))
cc_orig_dir   <- sapply(names(data_itn$itn_d), function(y) get_cc(data_itn$itn_d[[y]], max_norm = normalise_cc, directed = TRUE))

## (1) get the number of surrogates needed
N_s <- length(readRDS(file = list.files(results_dir, full.names = TRUE)[1])$cc)

## shuffle edge weights
cc_shuffle_undir <- t(sapply(names(data_itn$itn_u), function(y) calc_stats(replicate(N_s, get_cc(shuffle_weights(data_itn$itn_u[[y]]), max_norm = normalise_cc, directed = FALSE)), y)))
cc_shuffle_dir   <- t(sapply(names(data_itn$itn_d), function(y) calc_stats(replicate(N_s, get_cc(shuffle_weights(data_itn$itn_d[[y]]), max_norm = normalise_cc, directed = TRUE)), y)))

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

suffix <- ifelse(normalise_cc, "c", "k")

saveRDS(list("cc_orig_undir"     = cc_orig_undir,
             "cc_orig_dir"       = cc_orig_dir,
             "cc_shuffle_undir"  = cc_shuffle_undir,
             "cc_shuffle_dir"    = cc_shuffle_dir,
             "cc_undir_fixed"    = cc_undir_fixed,
             "cc_undir_interval" = cc_undir_interval,
             "cc_dir_fixed"      = cc_dir_fixed,
             "cc_dir_interval"   = cc_dir_interval),
        file = file.path(output_dir, paste0("itn_results_", suffix, ".rds")), compress = "xz")

## -----------------------------------------------------------------------------
