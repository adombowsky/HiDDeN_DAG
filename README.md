# HiDDeN
This is the GitHub repository for "Learning Discrete Bayesian Networks with Hierarchical Dirichlet Shrinkage". In this article, we propose a novel framework for Bayesian inference of discrete Bayesian networks (DBNs), called the Hierarchical Directed Dirichlet Network (HiDDeN). All results can be reproduced using the code in this repository. 

## Content
All simulation studies and applications are contained in the main folder of the repository. The folders ```r``` and ```rcpp``` contain functions written in R and Rcpp. There is some overlap between functions in these two folders, and we add the suffix ```_rcpp``` to Rcpp functions that have an equivalent in the ```r``` folder. The LUCAS and METABRIC datasets are publicly available. 

### Experiments
* Section 5.1: ```parameter_learning_simulation_study.R```.
* Section 5.2: ```lucas_blankets.R```.
* Section 5.3: ```evaluating_two_dags.R```.

### Application
* Section 6: ```metabric.R```.

### Supplement
* Section E: ```g1_results.R``` and ```g2_results.R```.
* Section F.2: ```alarm_probabilities.R```.
* Section F.3: ```alarm_network.R```.

## Key Functions
The following are key functions that implement some of the algorithms in the text.
* ```r/mala_within_gibbs.R/mala_within_gibbs_sampler```: Algorithm 1 with a fixed parent set.
* ```r/mala_within_gibbs.R/classifier_MCMC_samples```: Computes MCMC samples needed to evaluate equation (14) with a fixed parent set.
* ```rcpp/mala.cpp/mala_edge_sampler_rcpp```: Algorithm 1 + Algorithm 2; the MALA-within-Gibbs sampler where edges are represented via indicator functions.
* ```rcpp/mala.cpp/mala_model_sampler_rcpp```: Algorithm 1 + Algorithm 3; the MALA-within-Gibbs sampler that takes candidate parents/DAGs as input.
* ```rcpp/mala.cpp/log_posterior_samples```: Computes MCMC samples of the log-posterior.