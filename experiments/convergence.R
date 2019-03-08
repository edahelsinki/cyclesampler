## -----------------------------------------------------------------------------
## Perform the first convergence experiment
##
##       Sample networks where the vertex weights are preserved exactly,
##       while the edge weights w(e) are allowed to vary, subject to
##       interval constraints.
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla convergence.R 2000 100 1 datadir outputdir 0 1
##         
## -----------------------------------------------------------------------------

library(cyclesampler)
source("utilities.R")
source("init.R")

## -----------------------------------------------------------------------------
## Get command-line arguments
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

## -----------------------------------------------------------------------------
## Set the random seed
## -----------------------------------------------------------------------------
set.seed(42)

## -----------------------------------------------------------------------------
## Get the dataset
## -----------------------------------------------------------------------------
n_steps      <- as.numeric(args[1])                  ## number of total steps to take
stepsize     <- as.numeric(args[2])                  ## how often to save
ds           <- dataset_list[[as.numeric(args[3])]]  ## index of the dataset

datadir      <- args[4]                              ## data directory
outputdir    <- args[5]                              ## output directory

running_time <- as.numeric(args[6])                  ## estimate running time (Boolean)
convergence  <- as.numeric(args[7])                  ## estimate convergence (Boolean)
## -----------------------------------------------------------------------------

data <- read_data(name = ds$fname, basepath = datadir)

res  <- c(ds, get_data_properties(data))

out <- list()

out$t_normalize <- system.time( data <- normalizedata(data, scale_values = FALSE) )[3]
out$t_init      <- system.time( X    <- cyclesampler(data, a = rep(min(data[, 3]), nrow(data)), b = rep(max(data[, 3]), nrow(data)) ))[3]

## -----------------------------------------------------------------------------
## -- properties
## -----------------------------------------------------------------------------

res <- c(res, "n_nodes" = max(data[, 2]))
res <- c(res, "n_cycles" = X$getncycles())
saveRDS(res, file = paste0(outputdir, "/", ds$name, "_prop.rds"), compress = "xz")

## -----------------------------------------------------------------------------
## -- running time
## -----------------------------------------------------------------------------
##
if (running_time) {
    out$t_sample         <- replicate(10, system.time( tmp  <- X$samplecycles2(n = X$getncycles()) )[3])
    out$t_sample_average <- mean(out$t_sample)
    saveRDS(out, file = paste0(outputdir, "/", ds$name, "_running_time.rds"), compress = "xz")
}

## -----------------------------------------------------------------------------
## -- convergence
## save the dataset after having taken <stepsize> steps
## -----------------------------------------------------------------------------
if (convergence) {
    nullstate <- X$getstate()
    resvec    <- rep(NA, n_steps)

    attr(resvec, "n_steps")       <- n_steps
    attr(resvec, "stepsize")      <- stepsize
    attr(resvec, "ds")            <- ds
    attr(resvec, "time_start")    <- Sys.time()

    for (i in seq.int(n_steps / stepsize)) {
        a <- (i-1) * stepsize + 1
        b <- i * stepsize

        attr(resvec, "t0_current") <- Sys.time()

        resvec[a:b] <- estimate_convergence(X, s_orig = nullstate, n_samples = stepsize, ind = seq.int(nrow(data)))

        attr(resvec, "t1_current") <- Sys.time()
        attr(resvec, "time_stop")  <- Sys.time()

        saveRDS(resvec, file = paste0(outputdir, "/", ds$name, "_convergence.rds"), compress = "xz")
    }
}
## -----------------------------------------------------------------------------

