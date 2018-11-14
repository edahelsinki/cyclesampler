## -----------------------------------------------------------------------------
## Preprocess the datasets
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla preprocess_datasets.R path_to_dataset_directory
##         
## -----------------------------------------------------------------------------

source("utilities.R")

## --------------------------------------------------------------------------------
## Get command-line arguments
## --------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
bp   <- args[1]                               ## path to the dataset directory

## --------------------------------------------------------------------------------
## Preprocess last.fm
## --------------------------------------------------------------------------------
data <- read_data("lastfm_listening_count.csv", basepath = bp)
data <- preprocess_data(data, cutoff = 2500)
write.table(data, file = file.path(bp, "lastfm_listening_count_trimmed.csv"), sep = ";", row.names = FALSE, col.names = FALSE)

## --------------------------------------------------------------------------------
## Preprocess TasteProfile
## --------------------------------------------------------------------------------
data <- read_data("tasteprofile_processed.csv", basepath = bp)
data <- preprocess_data(data, cutoff = 20)
write.table(data, file = file.path(bp, "tasteprofile_processed_trimmed.csv"), sep = ";", row.names = FALSE, col.names = FALSE)

## --------------------------------------------------------------------------------
## Preprocess Book Crossing (BX-Books)
## --------------------------------------------------------------------------------
data <- read_data(name = "BX-Book-Ratings-processed.csv", basepath = bp)
data <- data[data[,3] > 0, ]              ## only use non-zero entries (explicit ratings)
ind  <- which(duplicated(data[, 1:2]))    ## remove duplicates
data <- data[-ind,]
write.table(data, file = file.path(bp, "BX-Book-Ratings-processed_trimmed.csv"), sep = ";", row.names = FALSE, col.names = FALSE)

## --------------------------------------------------------------------------------
