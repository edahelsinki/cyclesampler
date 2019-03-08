## --------------------------------------------------
## Retrieve and preprocess the
##    "Expanded Trade and GDP Data"
##    (International Trade Network dataset)
##
##    http://ksgleditsch.com/exptradegdp.html
## --------------------------------------------------
##
## Usage:
##
##        Rscript --vanilla get_itn_data.R destdir
##
## This will download the data, preprocess it
## and create an rds-file containing the trade network.
##
## The argument 'destdir' specifies the destination directory
## --------------------------------------------------

## --------------------------------------------------
## Command-line arguments
## --------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1)
    stop("Argument 'destdir' not provided!")

destdir <- args[1]


## --------------------------------------------------
## Libraries
## --------------------------------------------------

library(plyr)

## --------------------------------------------------
## Helper function
## --------------------------------------------------

## Map country identifiers to continuous range and convert to matrices
## perform normalisation using average edge weight
map_vals <- function(data, scalefactor = NULL) {
    data <- unname(as.matrix(data))

    ## remove zero-weight edges
    ind <- which(data[, 3] == 0)

    if (length(ind) > 0)
        data <- data[-ind, ]

    ## normalize by average weight
    data[, 3] <- data[, 3] / mean(data[, 3])

    ## remove year column
    data <- data[, -4]

    from_vertices <- sort(unique(c(data[, 1], data[, 2])))
    to_vertices   <- seq.int(length(from_vertices))

    data[, 1] <- plyr::mapvalues(x = data[, 1], from = from_vertices, to = to_vertices)
    data[, 2] <- plyr::mapvalues(x = data[, 2], from = from_vertices, to = to_vertices)

    data
}

## --------------------------------------------------
## Retrieve the data
## --------------------------------------------------
dataURL <- "http://ksgleditsch.com/data/exptradegdpv4.1.zip"
download.file(dataURL, destfile = file.path(destdir, basename(dataURL)))
unzip(file.path(destdir, basename(dataURL)), exdir = destdir)

## --------------------------------------------------
## Create networks
## --------------------------------------------------
data <- read.csv(file.path(destdir, "expdata.asc"), sep = " ", stringsAsFactors = FALSE)

## --------------------
## (1) undirected network of exports
## --------------------
v    <- 0.5 * apply(data[, c("expab", "expba", "impab", "impba")], 1, FUN = sum)
tmp1 <- cbind(data[, c("numa", "numb")], v, data[, "year"])

itn_u <- as.matrix(tmp1)
colnames(itn_u) <- NULL

## Split dataset per year and preprocess
itn_u_yearwise <- split(as.data.frame(itn_u), itn_u[, 4])
itn_u_yearwise <- lapply(itn_u_yearwise, function(df) map_vals(df))

## --------------------

## --------------------
## (2) directed network of exports
## --------------------
ew_ab <- 0.5 * apply(data[, c("expab", "impba")], 1, FUN = sum)
ew_ba <- 0.5 * apply(data[, c("expba", "impab")], 1, FUN = sum)

tmp1 <- cbind(data[, c("numa", "numb")], ew_ab, data[, "year"])
tmp2 <- cbind(data[, c("numb", "numa")], ew_ba, data[, "year"])

itn_d <- rbind(as.matrix(tmp1), as.matrix(tmp2))
colnames(itn_d) <- NULL

## Split dataset per year and preprocess
itn_d_yearwise <- split(as.data.frame(itn_d), itn_d[, 4])
itn_d_yearwise <- lapply(itn_d_yearwise, function(df) map_vals(df))

## --------------------------------------------------
## Save
## --------------------------------------------------

## Save networks
saveRDS(list("itn_u" = itn_u_yearwise, "itn_d" = itn_d_yearwise), file = "itn_data.rds", compress = "xz")

## --------------------------------------------------


