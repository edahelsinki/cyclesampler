Experiments
===========

The files in this folder are used to run the experiments for the paper
"Randomization Algorithms for Large Sparse Matrices"


Obtaining the data
-------------------

The directory 'data' contains bash-scripts that can be used to download the datasets.
Each of these scripts are run as, e.g., `./get_movielens.sh /tmp/movielens` where the only argument
is to some suitable temporary directory.

There are also Python (version 2) scripts and an R-script ("preprocess_datasets.R") for preprocessing
some of the datasets.


Examples and Experiments
-------------------------
The experiments are contained in the following files:

- **convergence.R** Sample networks where the vertex weights are preserved exactly, while the edge weights w(e) are allowed to vary, subject to interval constraints.
- **convergence_nwfix.R** Sample networks where both edge and vertex weights are allowed to vary, subject to interval constraints.

The examples are contained in the files **example_1.R** and **example_2.R**.


Analysis
---------
The convergence and runtime experiments are produced using the file **analyse_results.R**.
