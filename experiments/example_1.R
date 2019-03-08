## -----------------------------------------------------------------------------
## Example 1
##
##       Sample networks describing phone calls between individuals,
##       where edge and node weights are restricted to the
##       interval [0, 24] hours.
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla example_1.R
##
## -----------------------------------------------------------------------------

## --------------------------------------------------
## Load libraries and set the random seed
## --------------------------------------------------
library(cyclesampler)
source("utilities.R")
set.seed(42)

## --------------------------------------------------
## Helper functions
## --------------------------------------------------

## Helper function to get node weights This function does not take
## self-loops into account, by truncating the value-vector to the
## length of the original data matrix without self-loops.
gnw_helper <- function(data, values) {
    data[, 3] <- values[1:nrow(data)]
    get_node_weights(data)
}

## --------------------------------------------------
## Define the network
## --------------------------------------------------

data <- matrix(c(1, 2, 1.5,
                 1, 3, 5,
                 1, 6, 7,
                 2, 3, 4,
                 3, 4, 3,
                 4, 5, 8,
                 4, 6, 6),
               byrow = TRUE, ncol = 3)


n_edges    <- nrow(data)
n_vertices <- length(unique(c(data[, c(1, 2)])))

minval <- 0
maxval <- 24

tmp <- addselfloops(data, A = rep(minval, n_vertices), B = rep(maxval, n_vertices))

data_all <- rbind(data, tmp$data)

## Allowed edge weights in [0, 24] hours
a <- rep(minval, n_edges)
b <- rep(maxval, n_edges)

## range for self-loops
al <- c(a, tmp$a)
bl <- c(b, tmp$b)

## Sample N samples and calculate range of edge and node weights
N <- 10000

## ---------- CycleSampler ----------

## (1) Cyclesampler (node weights exact)
## (2) Cyclesampler (node and edge weights can vary on an interval

X1 <- cyclesampler(data, a = a, b = b)
X2 <- cyclesampler(data_all, a = al, b = bl)

run_sampler <- function(s) {
    s$samplecycles2(1000)
    s$getstate()
}

## edge weights
CC_ew_1 <- t(replicate(N, run_sampler(X1)))
CC_ew_2 <- t(replicate(N, run_sampler(X2)))

## node weights
CC_nw_1 <- apply(CC_ew_1, 1, function(i) gnw_helper(data = data, values = i))
CC_nw_2 <- apply(CC_ew_2, 1, function(i) gnw_helper(data = data, values = i))

## ---------- MaxEnt -----------------
data_me      <- data
data_me[, 3] <- data[, 3] / 24

Y   <- maxentsampler(data_me)
tmp <- Y$optimlambda(tol=0)

## Obtain 10000 samples from weight of vertex 2
ME_ew_1 <- t(replicate(N, Y$sample()) * 24)
ME_nw_1 <- apply(ME_ew_1, 1, function(i) gnw_helper(data = data, values = i))

## -----------------------------------

out <- matrix(nrow = 2, ncol = 8)
out[1, 1:2] <- range(CC_ew_1)
out[1, 3:4] <- range(CC_nw_1)

out[1, 5:6] <- range(CC_ew_2[, 1:7])
out[1, 7:8] <- range(CC_nw_2)

out[2, 1:2] <- range(ME_ew_1)
out[2, 3:4] <- range(ME_nw_1)

rownames(out) <- c("CycleSampler", "MaxEnt")

library(xtable)
print(xtable(out))

## --------------------------------------------------
