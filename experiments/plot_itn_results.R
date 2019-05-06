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
## -----------------------------------------------------------------------------
## Example usage:
##
##        Rscript --vanilla plot_itn_results.R <results_dir> <figure_dir>
##
## results_dir : path to the directory holding the results of the ITN analyses
## figure_dir  : path to the directory where the result figures should be saved
## -----------------------------------------------------------------------------

source("utilities.R")
library(latex2exp)
library(ggplot2)
library(scales)

## -----------------------------------------------------------------------------
## Get command-line arguments
## -----------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2)
    stop("Not enough arguments")

results_dir <- args[1]
figure_dir  <- args[2]

## -----------------------------------------------------------------------------
## Read the results
## -----------------------------------------------------------------------------
res_itn_c   <- readRDS(file.path(results_dir, "itn_results_c.rds"))
res_itn_k   <- readRDS(file.path(results_dir, "itn_results_k.rds"))

## -----------------------------------------------------------------------------
## Plot results
## -----------------------------------------------------------------------------

get_ylims <- function(res) {
    list("undir" = range(c(res$cc_shuffle_undir[, 2], res$cc_undir_fixed[, 3], res$cc_undir_interval[, 3])),
         "dir"   = range(c(res$cc_shuffle_dir[, 2], res$cc_dir_fixed[, 3], res$cc_dir_interval[, 3])))
}

plot_helper <- function(res, xlab, ylab, sci, main_undir, main_dir) {
    ylims <- get_ylims(res)

    list("p1" = make_itn_plot(res$cc_orig_undir, res$cc_shuffle_undir, res$cc_undir_fixed, res$cc_undir_interval, ylim = ylims$undir, main = main_undir, xlab = xlab, ylab = ylab, sci = sci, notitle = FALSE),
         "p2" = make_itn_plot(res$cc_orig_dir, res$cc_shuffle_dir, res$cc_dir_fixed, res$cc_dir_interval, ylim = ylims$dir, main = main_dir, xlab = xlab, ylab = ylab, sci = sci, notitle = FALSE))
}

plots_c <- plot_helper(res_itn_c, xlab = "year", ylab = NULL, sci = TRUE,  main_undir = TeX("$\\bar{C}$ for undirected networks"), main_dir = TeX("$\\bar{C}$ for directed networks"))
plots_k <- plot_helper(res_itn_k, xlab = "year", ylab = NULL, sci = FALSE, main_undir = TeX("$\\bar{K}$ for undirected networks"), main_dir = TeX("$\\bar{K}$ for directed networks"))

ggsave(plot = plots_c$p1, filename = file.path(figure_dir, "c_undir.pdf"), width = 150, height = 100, units = "mm")
ggsave(plot = plots_c$p2, filename = file.path(figure_dir, "c_dir.pdf"),   width = 150, height = 100, units = "mm")
ggsave(plot = plots_k$p1, filename = file.path(figure_dir, "k_undir.pdf"), width = 150, height = 100, units = "mm")
ggsave(plot = plots_k$p2, filename = file.path(figure_dir, "k_dir.pdf"),   width = 150, height = 100, units = "mm")

## -----------------------------------------------------------------------------
