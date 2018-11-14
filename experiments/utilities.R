## ======================================================================
## Utilities for data reading and for getting data properties
## ======================================================================

##' Read data
##'
##' @param fname Full path to data file (in csv format) to read.
##'
##' @return The data as a matrix, with duplicated items removed (in columns 1 and 2).
##'
##' @export
read_data <- function(name, basepath) {
    tmp <- as.matrix(read.csv(file.path(basepath, name), sep = ";", header = FALSE))
    colnames(tmp) <- NULL
    tmp
}


##' Get properties of the data.
##'
##' @param triplets Rating matrix
##' @param values Values to be used, replacing the third column in the matrix triplets.
##'
##' @return List of properties describing the data (rowsums, colsums, ...).
##'
##' @export
get_data_properties <- function(triplets, values = NULL) {
    if (! is.null(values))
        triplets[, 3] <- values

    nr <- length(unique(triplets[, 1]))
    nc <- length(unique(triplets[, 2]))

    list("nrow"    = nr,
         "ncol"    = nc,
         "n_edges" = nrow(triplets),
         "density" = exp(log(length(triplets[, 3])) - (log(nr) + log(nc))),
         "v_min"   = min(triplets[, 3]),
         "v_max"   = max(triplets[, 3]))
}

## ======================================================================
## Utilities for estimating the convergence
## ======================================================================

#' Calculate the Frobenius norm between two vectors
#'
#' @param v1 A vector
#' @param v2 A vector
#'
#' @return The Frobenius norm between v1 and v2.
#'
#' @export
frobenius_norm <- function(v1, v2) {
    sqrt(sum((v1 - v2)^2))
}


#' Run the sampler for n_steps and return the current state.
#'
#' @param X A cyclesampler object.
#' @param n_steps The number of steps to run the cycler.
#'
#' @return State of the sampler after running it.
#'
#' @export
sample_and_return_state <- function(X, n_steps) {
    X$samplecycles2(n_steps)
    X$getstate()
}


#' Return the development of the Frobenius norm when running the sampler.
#'
#' Calculates the Frobenius norm with respect to the state when calling the function.
#' The number of cycles is used as the thinning.
#'
#' @param X A cyclesampler object.
#' @param n_samples The number of samples to obtain
#' @param ind The indices of the true edges in the graph (used to exclude self-loops)
#'
#' @return The Frobenius norm for each of the n_samples with respect to the starting state.
#'
#' @export
estimate_convergence <- function(X, s_orig, n_samples = 1000, ind) {
    ## s_orig <- X$getstate()
    sapply(seq.int(n_samples), function(i) frobenius_norm(s_orig[ind], sample_and_return_state(X, X$getncycles() )[ind] ))
}


## ======================================================================
## Utilities for preprocessing
## ======================================================================

#' Remove all triplets <user, item, rating> where <rating> exceeds cutoff.
#'
#' @param data Matrix
#' @param cutoff Cutoff value for the rating.
#'
#' @return The data with some triplets removed.
#'
#' @export
preprocess_data <- function(data, cutoff) {
    ind_remove <- (data[, 3] > cutoff)
    data[-which(ind_remove), ]
}

## ======================================================================
## Utilities for plotting the results from the convergence experiments
## ======================================================================


#' Read experimental results and scale the Frobenius norm to
#' the interval [0, 1]
#'
#' @param basepath The location of the results (rds files)
#' @param fpattern Pattern to use for matching filenames.
#' @param log_downsample Should the data be downsampled using a
#'     logarithmic spacing, useful for visualisation purposes
#'     (Boolean, defailt is \code{FALSE}).
#' @param jitter Should the data be jittered to prevent overplotting
#'     (Boolean, defailt is \code{FALSE}).
#'
#' @return A list with the data as a matrix, as a data frame, a vector
#'     containing the current sample for each dataset and a list
#'     containing the indices of the data samples.
#'
#' @export
read_data_convergence <- function(basepath, fpattern, log_downsample = FALSE, jitter = FALSE) {
    flist <- list.files(path = basepath, pattern = fpattern, full.names = TRUE)

    out <- sapply(flist, function(f) readRDS(f), simplify = "array")
    colnames(out) <- gsub("_", " ", gsub(fpattern, "", sapply(colnames(out), basename, USE.NAMES = FALSE)))

    ## Current last sample for all datasets
    cs       <- colSums(! is.na(out))

    ## Normalize the convergence times
    out_n <- sweep(out, 2, out[1, ], FUN = '-')
    out_n <- sweep(out_n, 2, apply(out_n, 2, max, na.rm = TRUE), FUN = '/')

    if (jitter) {
        for (i in seq.int(-3, 3))
            out_n[, (i + 4)] <- out_n[, (i + 4)] + i/100
    }

    if (log_downsample) {
        ind <- unique(floor(pracma::logseq(1, nrow(out_n), 1000)))
        out_n <- out_n[ind, ]
    } else {
        ind <- seq.int(nrow(out_n))
    }


    out_df   <- as.data.frame(out_n)
    out_df$t <- ind
    out_df   <- melt(out_df, id = "t")

    out_df$variable <- ordered(out_df$variable, levels = c("Last.fm", "MovieLens 100k", "FineFoods", "MovieLens 1M", "BookCrossing", "MovieLens 20M", "TasteProfile"))

    list("out_n" = out_n, "out_df" = out_df, "cs" = cs, "ind" = ind)
}

#' Read experimental results (Frobenius norm) from the convergence
#' experiments and visualise them.
#'
#' @param basepath The location of the results (rds files)
#' @param fpattern Pattern to use for matching filenames.
#' @param log_ds Should the data be downsampled using a
#'     logarithmic spacing, useful for visualisation purposes
#'     (Boolean, defailt is \code{FALSE}).
#'
#' @return A ggplot2 object.
#'
#' @export
load_results_and_plot <- function(basepath, fpattern, log_ds = TRUE) {
    res    <- read_data_convergence(basepath, fpattern, log_downsample = log_ds, jitter = TRUE)
    out_n  <- res$out_n
    out_df <- res$out_df
    cs     <- res$cs

    ## Plot the results
    point <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

    p <- ggplot(as.data.frame(out_df))

    ## p <- p + geom_vline(xintercept = cs[["MovieLens 20M"]], colour = brewer.pal(n = 7, name = "Set2")[6], linetype = "dashed")
    ## p <- p + geom_vline(xintercept = cs[["TasteProfile"]], colour = brewer.pal(n = 7, name = "Set2")[7], linetype = "dashed")

    p <- p + geom_line(aes(x = t, y = value, group = variable, colour = variable), size = 0.7)
    p <- p + scale_x_log10(breaks = c(0, 100, 1000, 10000, 100000), labels = point, minor_breaks = NULL)
    p <- p + scale_y_continuous(breaks = seq.int(0, 1, 0.25), minor_breaks = NULL)

    p <- p + scale_colour_brewer(palette = "Set2")


    p <- p + ylab("Normalised Frobenius norm")
    p <- p + xlab("steps / null space dimensionality")
    p <- p + theme_bw()

    p <- p + theme(legend.position=c(0.59, 0.28),
                   legend.title = element_blank(),
                   legend.background = element_rect(fill = "white"),
                       ## element_blank(),
                   legend.key.height=unit(0.7, "line")
                   )
    p <- p + guides(colour = guide_legend(override.aes = list(size = 3)))

    p
}


#' Helper function to be used with xtable to bold row names.
#'
#' @param x String to bold.
#'
#' @return The string x wrapped in the LaTeX-command for bold.
#'
#' @export
bold <- function(x) {
    paste0('\\textbf{', x, '}')
}


#' Create a table with the properties of the networks (size, sampling
#' time etc).
#' 
#'
#' @param respath The location of the results (rds files)
#'
#' @return A table with the results.
#'
#' @export
make_network_property_table <- function(respath) {

    out <- NULL
    dsname <- c()

    for (ds in dataset_list) {
        print(ds$name)
        f1 <- file.path(respath, paste0(ds$name, "_prop.rds"))
        f2 <- file.path(respath, paste0(ds$name, "_running_time.rds"))
        f3 <- file.path(respath, paste0(ds$name, "_running_time_nwfix.rds"))

        if (file.exists(f1) & (file.exists(f2)) & (file.exists(f3))) {
            res1 <- readRDS(f1)
            res2 <- readRDS(f2)
            res3 <- readRDS(f3)

            dens <- exp( log(res1$n_edges) - (log(res1$nrow) + log(res1$ncol)))
            out <- rbind(out, c(res1$nrow,
                                res1$ncol,
                                dens,
                                res1$n_edges,
                                res1$n_nodes,
                                res1$n_cycles,

                                res2$t_normalize + res2$t_init,
                                res2$t_sample,

                                res3$t_normalize + res3$t_addselfloops + res3$t_init,
                                res3$t_sample))
            dsname <- c(dsname, ds$name)
        }

    }

    rownames(out) <- dsname
    colnames(out) <- c("nrow", "ncol", "density", "n_edges", "n_nodes", "n_cycles", "t_init_1", "t_sample_1", "t_init_2", "t_sample_2")

    ## order the data in increasing order of size
    out <- out[order(out[, colnames(out) == "n_edges"]), ]

    ## print the results as a latex table
    ## print(xtable(out), booktabs = TRUE)
    out
}
