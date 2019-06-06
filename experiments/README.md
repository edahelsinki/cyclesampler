Experiments
===========

The files in this folder are used to run the experiments for the following paper:

> Kai Puolamäki, Andreas Henelius, Antti Ukkonen. Randomization algorithms for large sparse networks.
> Physical Review E 99, 053311, 2019. <https://doi.org/10.1103/PhysRevE.99.053311>


Obtaining the data
-------------------

The directory 'data' contains bash-scripts and an R-script that can be used to download the datasets.
Each of the bash-scripts are run as, e.g., `./get_movielens.sh /tmp/movielens` where the only argument
is to some suitable temporary directory. The R-script is run as, e.g., `Rscript --vanilla get_itn_data.R /tmp/itndata`.

There are also Python (version 2) scripts and an R-script ("preprocess_datasets.R") for preprocessing
some of the datasets.


Examples and Experiments
-------------------------
The examples are contained in the files **example_1.R** and **example_2.R**.

The experiments on bipartite recommender datasets are contained in the following files:

- **convergence.R** Sample networks where the vertex strengths are preserved exactly, while the edge weights w(e) are allowed to vary subject to interval constraints.
- **convergence_nwfix.R** Sample networks where both edge weights and vertex strengths are allowed to vary subject to interval constraints.
- **analyse_results.R** Used to produce the convergence and runtime experiments, after first having run convergence.R and convergence_nwfix.R.

The experiments on the International Trade Network (ITN) data is contained in the following files:

- **analyse_itn.R** Calculate the average clustering coefficient for the ITN networks.
- **analyse_itn_results.R** Produce plots with confidence bands for the average clustering coefficient for surrogates for the ITN data, after first having run analyse_itn.R.

