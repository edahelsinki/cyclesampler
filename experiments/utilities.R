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

#' Calculate the L2 norm between two vectors
#'
#' @param v1 A vector
#' @param v2 A vector
#'
#' @return The L2 norm between v1 and v2.
#'
#' @export
l2_norm <- function(v1, v2) {
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


#' Return the development of the L2 norm when running the sampler.
#'
#' Calculates the L2 norm with respect to the state when calling the function.
#' The number of cycles is used as the thinning.
#'
#' @param X A cyclesampler object.
#' @param n_samples The number of samples to obtain
#' @param ind The indices of the true edges in the graph (used to exclude self-loops)
#'
#' @return The L2 norm for each of the n_samples with respect to the starting state.
#'
#' @export
estimate_convergence <- function(X, s_orig, n_samples = 1000, ind) {
    sapply(seq.int(n_samples), function(i) l2_norm(s_orig[ind], sample_and_return_state(X, X$getncycles() )[ind] ))
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
## Utilities for the analysis of the International Trade Networks
## ======================================================================

##' Get a sampler based on the data.
##'
##' @param data A n x 3 matrix containing triplets representing the
##'     network
##' @param normalize Normalize the data. Boolean. Values are never scaled when normalizing.
##' @param sampler_type A string denoting the type of constraints to use
##'     on vertex strengths. Either "fixed" (preserve strengths
##'     exactly) or "interval" (preserve strengths on an interval).
##' @param interval The interval for preserving vertex strengths,
##'     if sampler_type is interval. Default is interval = c(0.90,
##'     1.10), i.e., [0.9 * ## nw, 1.1 * nw], allowing a +-10 percent
##'     variation in node weights.
##'
##' @return The data as a matrix, with duplicated items removed (in columns 1 and 2).
##'
##' @export
get_sampler <- function(data, normalize, sampler_type, interval = c(0.90, 1.10)) {
    ## Normalize the range of data values
    if (normalize) {
        data_n <- normalizedata(data, scale_values = FALSE)
    } else {
        data_n <- data
    }

    ## Preserve vertex strengths exactly
    if (sampler_type == "fixed") {
        ew      <- matrix(nrow = nrow(data_n), ncol = 2)
        ew[, 1] <- min(data_n[, 3])
        ew[, 2] <- max(data_n[, 3])
    }

    ## Preserve vertex strengths on an interval
    if (sampler_type == "interval") {
        nw <- get_node_weights(data_n)

        sl_tmp <- addselfloops(data_n, A = (nw * interval[1]), B = (nw * interval[2]))

        ew      <- matrix(nrow = (nrow(data_n) + length(sl_tmp$a)), ncol = 2)
        ew[, 1] <- c(rep(min(data_n[, 3]), nrow(data_n)), sl_tmp$a)
        ew[, 2] <- c(rep(max(data_n[, 3]), nrow(data_n)), sl_tmp$b)

        data_n <- rbind(data_n, sl_tmp$data)
    }

    ## Return (i) sampler,  (ii) indices of original (non-self-loop) data items and (iii) the network
    list("sampler"      = cyclesampler(data_n, a = ew[, 1], b = ew[, 2]),
         "ind"          = seq.int(nrow(data)),
         "data"         = data,
         "data_n"       = data_n,
         "ew"           = ew,
         "normalize"    = normalize,
         "sampler_type" = sampler_type)
}


##' Calculate average clustering coefficient for a matrix.
##'
##' @param data A n x 3 matrix containing triplets representing the
##'     network
##' @param max_norm Scale the average clustering coefficient using the
##'     maximum edge weight. Boolean. Default is TRUE.
##' @param directed Boolean. Indicates if the network is directed or
##'     not.
##'
##' @return The average clustering coefficient.
##'
##' @export
get_cc <- function(data, max_norm = TRUE, directed = FALSE) {
    library(Matrix)
    n <- length(unique(c(data[, 1], data[, 2])))

    if (directed) {
        mat <- sparseMatrix(i = data[, 1], j = data[, 2], x = data[, 3])
    } else {
        mat <- sparseMatrix(i = data[, 1], j = data[, 2], x = data[, 3], symmetric = TRUE)
    }

    mat <- mat^(1/3)
    tmp <- (mat %*% mat %*% mat)
    cc  <- diag(tmp) / ((n - 1) * (n - 2))

    if (max_norm)
        cc <- cc / max(data[, 3])

    sum(cc) / sum(cc > 0)
}


##' Helper function to get clustering coefficient from a data sample.
##'
##' @param data A n x 3 matrix containing triplets representing the
##'     network.
##' @param state Current state of the sampler. Vector.
##' @param ind The indices of the state vector representing the
##'     non-self-loop edges in the network.
##' @param directed Boolean. Indicates if the network is directed or
##'     not.
##' @param max_norm Scale the average clustering coefficient using the
##'     maximum edge weight. Boolean. Default is TRUE.
##'
##' @return The average clustering coefficient.
##'
##' @export
get_cc_helper <- function(data, state, ind, directed, max_norm = TRUE) {
    get_cc(cbind(data[, c(1, 2)], state[ind]), directed = directed, max_norm = max_norm)
}


##' Helper function to get clustering coefficient from a sampler.
##'
##' @param X A CycleSampler.
##' @param data A n x 3 matrix containing triplets representing the
##'     network.
##' @param n_samples Number of sampler of the clustering coefficient to obtain.
##' @param n_thin Thinning used when sampling the chain.
##' @param ind The indices of the state vector representing the
##'     non-self-loop edges in the network.
##' @param directed Boolean. Indicates if the network is directed or
##'     not.
##' @param max_norm Scale the average clustering coefficient using the
##'     maximum edge weight. Boolean. Default is TRUE.
##'
##' @return The average clustering coefficient.
##'
##' @export
sample_cc <- function(X, data, n_samples, n_thin, ind, directed, max_norm = TRUE) {
    replicate(n_samples,
              get_cc_helper(data,
                            state = sample_and_return_state(X, n_steps = (n_thin * X$getncycles())),
                            ind = ind,
                            directed = directed,
                            max_norm = max_norm))
}


##' Sample n_samples using the Besag-Clifford 'serial' method
##' It is here assumed that the sampler has not been "burned-in".
##'
##' @param sampler The output from \code{get_sampler}.
##' @param n_samples Number of sampler of the clustering coefficient
##'     to obtain.
##' @param n_thin Thinning used when sampling the chain.
##' @param directed Boolean. Indicates if the network is directed or
##'     not.
##' @param max_norm Scale the average clustering coefficient using the
##'     maximum edge weight. Boolean. Default is TRUE.
##'
##' @return The average clustering coefficient.
##'
##' @export
get_cc_besag_clifford_serial <- function(sampler, n_samples, n_thin, directed, max_norm = TRUE) {
    ## get the sampler and store the starting state
    X  <- sampler$sampler
    x0 <- X$getstate()

    ## number of samples to obtain from the forward / backward chains
    ind <- sample(seq.int(n_samples), 1) 
    n_f <- ind - 1
    n_b <- n_samples - ind
    
    ## sample the forward chain chain at intervals of n_thin starting from x0
    cc_f <- sample_cc(X, data = sampler$data, n_samples = n_f, n_thin = n_thin, directed = directed, ind = sampler$ind, max_norm = max_norm)

    ## sample the backward chain chain at intervals of n_thin starting from x0
    X$setstate(x0)
    cc_b <- sample_cc(X, data = sampler$data, n_samples = n_b, n_thin = n_thin, directed = directed, ind = sampler$ind, max_norm = max_norm)

    ## return clustering coefficients
    c(cc_f, cc_b)
}


##' Randomly shuffle edge weights in a network.
##'
##' @param data A n x 3 matrix containing triplets representing the
##'     network.
##' @return The network with the edge weights randomly shuffled.
##'
##' @export
shuffle_weights <- function(data) {
    data[, 3] <- sample(data[, 3])
    data
}


##' Sample n_samples using the Besag-Clifford 'serial' method
##' It is here assumed that the sampler has not been "burned-in".
##'
##' @param sampler The output from \code{get_sampler}.
##' @param n_samples Number of sampler of the clustering coefficient to obtain.
##' @param n_thin Thinning used when sampling the chain.
##' @param directed Boolean. Indicates if the network is directed or
##'     not.
##'
##' @return The average clustering coefficient.
##'
##' @export
calc_stats <- function(x, year) {
    m <- mean(x)
    s <- sd(x)
    c(as.numeric(year), m - s, m + s)
}


##' Helper function to read results from the ITN analysis and return a
##' confidence interval of width 1 standard deviation.
##'
##' @param fname The RDS-file to read.
##' 
##' @return Return the mean and confidence intervals of width 1
##'     standard deviation.
##'
##' @export
read_itn_res_helper <- function(fname) {
    year <- as.numeric(strsplit(basename(fname), "_")[[1]][1])
    tmp  <- readRDS(fname)

    ## calculate mean and standard deviation
    x  <- tmp$cc
    m  <- mean(x)
    sd <- sd(x)
    c(year, m-sd, m+sd)
}
    

##' Read results from the ITN analysis and return the confidence
##' intervals of width 1 standard deviation for each year in the ITN
##' results.
##'
##' @param basepath The path to the directory holding the RDS-files.
##' @param network_type A string denoting the type of the network
##'     ("dir" or "undir").
##' @param sampler_type A string denoting the type of the sampler
##'     ("fixed" or "interval").
##' 
##' @return Return a matrix with confidence intervals of width 1
##'     standard deviation for each year in the ITn results.
##'
##' @export
read_itn_results <- function(basepath, network_type, sampler_type) {
    flist         <- list.files(path = basepath, pattern = paste0("_", network_type, "_", sampler_type), full.names = TRUE)
    res           <- t(sapply(flist, function(fname) read_itn_res_helper(fname)))
    rownames(res) <- NULL
    colnames(res) <- c("year", "ci_lower", "ci_upper")
    res
}


## ======================================================================
## Utilities for plotting the results from the convergence experiments
## ======================================================================

#' Read experimental results (L2 norm) from the convergence
#' experiments and visualise them.
#'
#' @param basepath The location of the results (rds files)
#' @param dsname Name of dataset to plot.
#' @param dsname Pattern to use for matching filenames.
#' @param log_ds Should the data be downsampled using a
#'     logarithmic spacing, useful for visualisation purposes
#'     (Boolean, defailt is \code{FALSE}).
#'
#' @return A ggplot2 object.
#'
#' @export
load_results_and_plot_l2_norm <- function(basepath, dsname, suffix, log_ds = TRUE) {
    ## read data
    res <- readRDS(file.path(basepath, paste0(dsname, "_", suffix,".rds")))

    ## downsample and normalise
    if (log_ds) {
        ind <- unique(floor(pracma::logseq(1, length(res), 1000)))
        res <- res[ind]
    } else {
        ind <- seq.int(length(res))
    }
        
    res <- res - res[1]
    res <- res / max(res, na.rm = TRUE)
    
    res_df <- data.frame("t" = ind, "value" = res)
    
    ## Plot the results
    point <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

    p <- ggplot(as.data.frame(res_df))

    p <- p + geom_line(aes(x = t, y = value), color = "black", size = 0.7)
    p <- p + scale_x_log10(breaks = c(1, 100, 1000, 10000, 100000), labels = point, minor_breaks = NULL) ##, expand = c(0.2, 0, 0.2, 0))
    p <- p + scale_y_continuous(breaks = seq.int(0, 1, 0.25), minor_breaks = NULL)

    p <- p + ylab(TeX("Normalised $l^2$ norm"))
    p <- p + xlab("cycle steps")
    p <- p + theme_bw()

    p <- p + geom_rect(aes(xmin = 1000+50, xmax = 100000-5000, ymin = 0.25+0.05, ymax = 0.50-0.05), fill = "white")
    p <- p + geom_text(aes(x = 100000-5000, y = 0.38), label = gsub("_", " ", dsname), hjust = "right")
    
    p <- p + theme(legend.position = "none")
    p <- p + theme(panel.border = element_blank())
    p <- p + theme(axis.line = element_line(colour = "black"))

    p <- p + theme(axis.text = element_text(size = 10, colour = "black"),
                   axis.title = element_text(size = 10, colour = "black"))

    
    p <- p + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
    
    p <- p + theme(plot.margin = unit(c(5,7,5,7), unit = "mm"))  # t r b l
    
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
                                res2$t_sample_average,
                                
                                res3$t_normalize + res3$t_addselfloops + res3$t_init,
                                res3$t_sample_average))
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

## ======================================================================
## Utilities for plotting the results from the analysis of the
## International Trade Networks.
## ======================================================================

#' Plot the ITN results.
#'
#' @param cc_orig Clustering coefficients from the original data
#' @param cc_shuffle Clustering coefficients from data with edge weights permuted at random.
#' @param cc_sampler1 Clustering coefficients from data with vertex strengths preserved exactly.
#' @param cc_sampler2 Clustering coefficients from data with vertex strengths preserved on an interval.
#' @param main Title of the plot.
#' @param ylims Y-axis limits.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param sci Should scientific notation be used on the y-axis?
#' @param notitle Boolean indicating íf the title should be plotted or not. Default is FALSE.
#'
#' @return A plot.
#'
#' @export
make_itn_plot <- function(cc_orig, cc_shuffle, cc_sampler1, cc_sampler2, main = NULL, ylims, xlab, ylab, sci = TRUE, notitle = FALSE) {
    df_orig     <- data.frame(year = as.numeric(names(cc_orig)), value = cc_orig)
    df_shuffle  <- data.frame(year = cc_shuffle[, 1], ci_lower = cc_shuffle[, 2], ci_upper = cc_shuffle[, 3])
    df_sampler1 <- data.frame(year = cc_sampler1[, 1], ci_lower = cc_sampler1[, 2], ci_upper = cc_sampler1[, 3])
    df_sampler2 <- data.frame(year = cc_sampler2[, 1], ci_lower = cc_sampler2[, 2], ci_upper = cc_sampler2[, 3])

    p <- ggplot()

    a <- 0.5

    p <- p + geom_line(data = df_shuffle, aes(x = year, y = ci_lower), colour = "black", linetype = "dotdash")
    p <- p + geom_line(data = df_shuffle, aes(x = year, y = ci_upper), colour = "black", linetype = "dotdash")

    p <- p + geom_line(data = df_sampler1, aes(x = year, y = ci_lower), colour = "black", linetype = "dashed")
    p <- p + geom_line(data = df_sampler1, aes(x = year, y = ci_upper), colour = "black", linetype = "dashed")

    p <- p + geom_line(data = df_sampler2, aes(x = year, y = ci_lower), colour = "black", linetype = "dotted")
    p <- p + geom_line(data = df_sampler2, aes(x = year, y = ci_upper), colour = "black", linetype = "dotted")
    
    p <- p + geom_line(data = df_orig, aes(x = year, y = value), size = 1)

    if (sci)
        p <- p + scale_y_log10(limits = ylims, breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
    else
        p <- p + scale_y_log10(limits = ylims, breaks = trans_breaks("log10", function(x) 10^x), labels = number_format(accuracy = 0.01))
    
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)

    p <- p + theme_bw()

    p <- p + theme(axis.text = element_text(size = 15),
                   axis.title = element_text(size = 20),
                   plot.title = element_text(hjust = 0.5, size = 20))

    p <- p + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

    p <- p + theme(panel.border = element_blank())
    p <- p + theme(axis.line = element_line(colour = "black"))

    if (! notitle)
        p <- p + ggtitle(main)
    
    p
}

## ======================================================================
