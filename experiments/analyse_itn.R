## -----------------------------------------------------------------------------
## 
## Analysis of the International Trade Network data
##
## Sample N_samples with a thinning of N_thin and calculate the clustering coefficient,
## using the Besag-Clifford 'serial' method
##
##       Different samplers are used to preserve (i) vertex strenghts exactly
##       and to presrve (ii) vertex strenghts on an interval
##
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla analyse_itn.R <data_type> <sampler_type> <year> <N_samples> <N_thin> <state_dir>
##
## -----------------------------------------------------------------------------

library(cyclesampler)
library(methods)
source("utilities.R")

## -----------------------------------------------------------------------------
## Get command-line arguments
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

data_type    <- args[1]                    ## "undir" or "dir"
sampler_type <- args[2]                    ## "fixed" or "interval

year         <- as.character(args[3])      ## 1948 to 2000

N_samples    <- as.numeric(args[4])        ## 1000
N_thin       <- as.numeric(args[5])        ## 1000

data_dir     <- as.character(args[6])      ## path to directory holding the datafile "itn_data.rds"
output_dir   <- as.character(args[7])      ## path to the directory where the results should be saved

## -----------------------------------------------------------------------------
## Set the random seed
## -----------------------------------------------------------------------------
set.seed(42)

## -----------------------------------------------------------------------------
## Read the data
## -----------------------------------------------------------------------------

data_itn <- readRDS(file.path(data_dir, "itn_data.rds"))

if (data_type == "undir") {
    data <- data_itn$itn_u[[year]]
    normalize <- FALSE
    directed  <- FALSE
}

if (data_type == "dir") {
    data <- data_itn$itn_d[[year]]
    normalize <- TRUE
    directed  <- TRUE
}

## -----------------------------------------------------------------------------
## Initialise sampler
## -----------------------------------------------------------------------------
sampler <- get_sampler(data = data, normalize = normalize, sampler_type = sampler_type, interval = c(0.75, 1.25))

## -----------------------------------------------------------------------------
## Acquire samples and calculate the clustering coefficient
## -----------------------------------------------------------------------------
t_start <- Sys.time()
res_cc  <- get_cc_besag_clifford_serial(sampler, n_samples = N_samples, n_thin = N_thin, directed = directed)
t_stop  <- Sys.time()

## -----------------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------------
saveRDS(list("cc"= res_cc,
             "cc_orig" = get_cc(data = data, max_norm = TRUE, directed = directed),
             "t_start" = t_start,
             "t_tstop" = t_stop),
        file = file.path(output_dir, paste0(year, "_", data_type, "_", sampler_type, ".rds")),
        compress = "xz")

## -----------------------------------------------------------------------------
