#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat expand_grid_cpp(List v) {
  Function eg("expand.grid");
  NumericMatrix x = internal::convert_using_rfunction(eg(v), "as.matrix");
  return(arma::mat(x.begin(), x.nrow(), x.ncol(), false));
}

// [[Rcpp::export]]
int sample_arma(arma::vec probs) {
  int K = probs.n_rows;
  IntegerVector clusts = Range(1,K);
  IntegerVector samp = RcppArmadillo::sample(clusts, 1, TRUE, probs);
  int s = samp(0);
  return(s);
}

// [[Rcpp::export]]
Rcpp::List concatenate_lists(Rcpp::List list1, Rcpp::List list2) {
  // create a new list with a size large enough to hold all elements
  Rcpp::List combined_list(list1.size() + list2.size());
  
  for (int i = 0; i < list1.size(); ++i) {
    combined_list[i] = list1[i];
  }
  
  for (int i = 0; i < list2.size(); ++i) {
    combined_list[i + list1.size()] = list2[i];
  }
  
  return combined_list;
}

// [[Rcpp::export]]
arma::vec zero_to_n(int n) {
  arma::vec x = arma::zeros(n+1,1);
  for (int i = 0; i<=n; i++) {
    x(i) = i;
  }
  return(x);
} 

// [[Rcpp::export]]
Rcpp::List list_of_categories(const arma::vec& x) {
  int p = x.n_elem;
  Rcpp::List category_list(p);
  for (int j=0; j<p; j++) {
    category_list[j] = zero_to_n(x(j));
  }
  return(category_list);
} 

// [[Rcpp::export]]
int cell_count(arma::mat x, arma::vec m) {
  int n = x.n_rows;
  int p = x.n_cols;
  int c = 0;
  for (int i=0; i<n;i++) {
    if (sum(x.row(i).t()==m)==p) {
      c=c+1;
    }
  }
  return(c);
}

// [[Rcpp::export]]
arma::vec vectorized_contingency(arma::mat x, arma::mat g) {
  // x = data
  // g = grid
  int G = g.n_rows;
  arma::vec ct = arma::zeros(G,1);
  for (int j=0; j<G; j++) {
    ct(j) = cell_count(x, g.row(j).t());
  }
  return(ct);
}

// [[Rcpp::export]]
arma::uvec remove_jth_element(arma::uvec v, int j) {
  // check if j is a valid 1-based index
  if (j < 0 || j > v.n_elem-1) {
    Rcpp::stop("Invalid index j. Must be between 0 and the vector length.");
  }
  // remove j
  arma::uvec v_prime(v.n_elem-1);
  if (j==0) {
    v_prime = v.subvec(j + 1, v.n_elem - 1);
  } else if (j==v.n_elem-1) {
    v_prime = v.subvec(0, j - 1);
  } else {
    v_prime = v.subvec(0, j - 1);
    v_prime = arma::join_cols(v_prime, v.subvec(j + 1, v.n_elem - 1)); 
  }
  return v_prime;
}

// [[Rcpp::export]]
arma::mat remove_jth_column(const arma::mat& X, int j) {
  // create a local, non-const copy of the matrix to modify
  arma::mat new_X = X;
  
  if (j >= new_X.n_cols) {
    Rcpp::stop("Error: Invalid column index. Index is out of bounds.");
  }
  
  new_X.shed_col(j);
  
  return new_X;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List create_grids_and_table_child_rcpp(const arma::vec& x_j, const arma::mat& x_Pa_j) {
  // number of categories
  int k_j = x_j.max();
  arma::rowvec k_Pa_j = arma::max(x_Pa_j, 0);

  // initialize grids 
  arma::mat par_grid;
  arma::mat par_child_grid;

  // create parent and parent-child grids
  par_grid = expand_grid_cpp(list_of_categories(k_Pa_j.t()));
  // create a temporary arma::vec from k_j to avoid error
  arma::vec k_j_vec(1);
  k_j_vec(0) = k_j;
  par_child_grid = expand_grid_cpp(list_of_categories(arma::join_cols(k_j_vec,k_Pa_j.t())));

  // create tables
  arma::vec mct = vectorized_contingency(x_Pa_j,par_grid);
  arma::vec ct = vectorized_contingency(join_rows(x_j,x_Pa_j),par_child_grid);

  // returning
  return Rcpp::List::create(Rcpp::Named("par_child_grid") = par_child_grid,
                            Rcpp::Named("par_grid") = par_grid,
                            Rcpp::Named("ct") = ct,
                            Rcpp::Named("mct") = mct);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List create_grids_and_table_root_rcpp(const arma::vec& x_j) {
  // number of categories
  int k_j = x_j.max();
  int n = x_j.n_rows;
  
  // create objects
  arma::mat par_grid = arma::zeros(1,1);
  arma::vec ks(2);
  ks(0) = k_j;
  ks(1) = 0;
  arma::mat par_child_grid = expand_grid_cpp(list_of_categories(ks));
  arma::vec mct(1);
  mct(0) = n;
  arma::vec ct = vectorized_contingency(join_rows(x_j,arma::zeros(n,1)),par_child_grid);
  
  // returning
  return Rcpp::List::create(Rcpp::Named("par_child_grid") = par_child_grid,
                            Rcpp::Named("par_grid") = par_grid,
                            Rcpp::Named("ct") = ct,
                            Rcpp::Named("mct") = mct);
}

// [[Rcpp::export]]
double log_node_score_rcpp(const arma::vec& t,
                           const arma::vec& ct,
                           const arma::vec& mct,
                           const arma::mat& par_grid,
                           const arma::mat& par_child_grid) {
  
  double beta_j = arma::sum(t);
  double A = 0.0;
  
  // normalizing constant
  for (arma::uword x = 0; x < par_grid.n_rows; ++x) {
    A += lgamma(beta_j) - lgamma(mct(x) + beta_j);
  }
  
  double B = 0.0;
  
  // node-parent kernel
  for (arma::uword x = 0; x < par_child_grid.n_rows; ++x) {
    B += lgamma(ct(x) + t(par_child_grid(x, 0))) - lgamma(t(par_child_grid(x, 0)));
  }
  
  return A + B;
}

// [[Rcpp::export]]
int update_edge_rcpp(const arma::vec& t,
                     const arma::uvec& z_minus,
                     const arma::vec& x_j,
                     const arma::vec& x_j_prime,
                     const arma::mat& X_minus,
                     double c_j,
                     double d_j) {
  
  arma::uword M = arma::sum(z_minus);
  arma::uword W = arma::sum(1 - z_minus);
  
  arma::mat x_Pa_j_minus;
  if (M > 0) {
    arma::uvec parent_indices = arma::find(z_minus == 1);
    x_Pa_j_minus = X_minus.cols(parent_indices);
  }
  
  // calculate the score for no-edge
  Rcpp::List obj_no_edge;
  if (M == 0) {
    obj_no_edge = create_grids_and_table_root_rcpp(x_j);
  } else {
    obj_no_edge = create_grids_and_table_child_rcpp(x_j, x_Pa_j_minus);
  }
  
  double w = log(d_j + W) + log_node_score_rcpp(
    t,
    Rcpp::as<arma::vec>(obj_no_edge["ct"]),
    Rcpp::as<arma::vec>(obj_no_edge["mct"]),
    Rcpp::as<arma::mat>(obj_no_edge["par_grid"]),
    Rcpp::as<arma::mat>(obj_no_edge["par_child_grid"])
  );
  
  // calculate the score for an edge
  Rcpp::List obj_with_edge;
  if (M == 0) {
    obj_with_edge = create_grids_and_table_child_rcpp(x_j, x_j_prime);
  } else {
    arma::mat joined_parents = arma::join_cols(x_Pa_j_minus.t(), x_j_prime.t()).t();
    obj_with_edge = create_grids_and_table_child_rcpp(x_j, joined_parents);
  }
  
  double m = log(c_j + M) + log_node_score_rcpp(
    t,
    Rcpp::as<arma::vec>(obj_with_edge["ct"]),
    Rcpp::as<arma::vec>(obj_with_edge["mct"]),
    Rcpp::as<arma::mat>(obj_with_edge["par_grid"]),
    Rcpp::as<arma::mat>(obj_with_edge["par_child_grid"])
  );
  
  // Sample z_j_jprime
  arma::vec log_probs(2);
  log_probs(0) = w;
  log_probs(1) = m;
  
  double max_log_prob = arma::max(log_probs);
  arma::vec probs = arma::exp(log_probs - max_log_prob);
  probs = probs / arma::sum(probs);
  
  Rcpp::IntegerVector values = Rcpp::IntegerVector::create(0, 1);
  return Rcpp::as<int>(Rcpp::sample(values, 1, false, Rcpp::NumericVector(probs.begin(), probs.end())));
}

// [[Rcpp::export]]
int update_edge_one_candidate_rcpp(const arma::vec& t,
                     const arma::vec& x_j,
                     const arma::vec& x_j_prime,
                     double c_j,
                     double d_j) {
  
  arma::uword M = 0;
  arma::uword W = 0;
  
  // calculate the score for no-edge
  Rcpp::List obj_no_edge = create_grids_and_table_root_rcpp(x_j);
  
  double w = log(d_j + W) + log_node_score_rcpp(
    t,
    Rcpp::as<arma::vec>(obj_no_edge["ct"]),
    Rcpp::as<arma::vec>(obj_no_edge["mct"]),
    Rcpp::as<arma::mat>(obj_no_edge["par_grid"]),
    Rcpp::as<arma::mat>(obj_no_edge["par_child_grid"])
  );
  
  // calculate the score for an edge
  Rcpp::List obj_with_edge = create_grids_and_table_child_rcpp(x_j, x_j_prime);
  
  double m = log(c_j + M) + log_node_score_rcpp(
    t,
    Rcpp::as<arma::vec>(obj_with_edge["ct"]),
    Rcpp::as<arma::vec>(obj_with_edge["mct"]),
    Rcpp::as<arma::mat>(obj_with_edge["par_grid"]),
    Rcpp::as<arma::mat>(obj_with_edge["par_child_grid"])
  );
  
  // Sample z_j_jprime
  arma::vec log_probs(2);
  log_probs(0) = w;
  log_probs(1) = m;
  
  double max_log_prob = arma::max(log_probs);
  arma::vec probs = arma::exp(log_probs - max_log_prob);
  probs = probs / arma::sum(probs);
  
  Rcpp::IntegerVector values = Rcpp::IntegerVector::create(0, 1);
  return Rcpp::as<int>(Rcpp::sample(values, 1, false, Rcpp::NumericVector(probs.begin(), probs.end())));
}


// [[Rcpp::export]]
arma::uvec update_parents_rcpp(const arma::vec& t,
                    arma::uvec& z,
                     const arma::vec& x_j,
                     const arma::mat& X_minus_j,
                     double c_j,
                     double d_j) {
  int M = z.n_elem;
  if (M==1) {
    z(0) = update_edge_one_candidate_rcpp(t,x_j,X_minus_j,c_j,d_j);
    } 
  else{
    for (int j_prime=0; j_prime<M; j_prime++) {
      z(j_prime) =
        update_edge_rcpp(t,
                         remove_jth_element(z,j_prime),
                         x_j,
                         X_minus_j.col(j_prime),
                         remove_jth_column(X_minus_j,j_prime),
                         c_j,
                         d_j);
    }
  }
  return(z);
}

// [[Rcpp::export]]
double log_h_rcpp(double t, double b, double rho, int k,
                  const arma::vec& N, const arma::vec& U) {
  
  double g = rho / (static_cast<double>(k) + 1.0);
  double result = -b * t + (g - 1.0) * log(t);
  double sum_lgamma_and_log_U = 0.0;
  for (arma::uword i = 0; i < N.n_elem; ++i) {
    sum_lgamma_and_log_U += lgamma(N(i) + t) - lgamma(t) + t * log(U(i));
  }
  return result + sum_lgamma_and_log_U;
}

// [[Rcpp::export]]
double grad_log_h_rcpp(double t, double b, double rho, int k,
                       const arma::vec& N, const arma::vec& U) {
  
  double g = rho / (static_cast<double>(k) + 1.0);
  double result = -b + (g - 1.0) / t;
  double sum_digamma_and_log_U = 0.0;
  for (arma::uword i = 0; i < N.n_elem; ++i) {
    sum_digamma_and_log_U += R::digamma(N(i) + t) - R::digamma(t) + log(U(i));
  }
  return result + sum_digamma_and_log_U;
}

// [[Rcpp::export]]
double mala_proposal_rcpp(double t_current, double step_size,
                          double b, double rho, int k,
                          const arma::vec& N, const arma::vec& U) {
  
  double gradient = grad_log_h_rcpp(t_current, b, rho, k, N, U);
  double t_proposal = t_current + (pow(step_size, 2) / 2.0) * gradient +
    step_size * arma::randn();
  return t_proposal;
}

// [[Rcpp::export]]
double log_mala_density_rcpp(double t_current, double t_proposal, double step_size,
                             double b, double rho, int k,
                             const arma::vec& N, const arma::vec& U) {
  
  double mean_val = t_current + (pow(step_size, 2) / 2.0) * grad_log_h_rcpp(t_current, b, rho, k, N, U);
  
  // Calculate the PDF and then take its natural logarithm
  return arma::log_normpdf(t_proposal, mean_val, step_size);
}

// [[Rcpp::export]]
Rcpp::List mala_update_rcpp(double t_current, double step_size,
                            double b, double rho, int k,
                            const arma::vec& N, const arma::vec& U) {
  
  // compute proposal
  double t_proposal = mala_proposal_rcpp(t_current, step_size, b, rho, k, N, U);
  
  // If proposal is negative, reject immediately
  if (t_proposal < 0) {
    return Rcpp::List::create(Rcpp::Named("next_state") = t_current,
                              Rcpp::Named("acceptance") = 0);
  } else {
    // calculate log acceptance ratio
    double log_acceptance_ratio =
      log_h_rcpp(t_proposal, b, rho, k, N, U) - log_h_rcpp(t_current, b, rho, k, N, U) +
      log_mala_density_rcpp(t_proposal, t_current, step_size, b, rho, k, N, U) -
      log_mala_density_rcpp(t_current, t_proposal, step_size, b, rho, k, N, U);
    
    // acceptance probability
    if (log(arma::randu()) < log_acceptance_ratio) {
      // Accept
      return Rcpp::List::create(Rcpp::Named("next_state") = t_proposal,
                                Rcpp::Named("acceptance") = 1);
    } else {
      // Reject
      return Rcpp::List::create(Rcpp::Named("next_state") = t_current,
                                Rcpp::Named("acceptance") = 0);
    }
  }
}

// [[Rcpp::export]]
double beta_update_rcpp(double beta, double n) {
  return R::rbeta(beta, n);
}

// [[Rcpp::export]]
Rcpp::List mala_within_gibbs_graph_rcpp(arma::vec t,
                                        const arma::vec& ct, const arma::vec& mct,
                                        const arma::mat& par_grid, const arma::mat& par_child_grid,
                                        double b, double rho, int k,
                                        const arma::vec& step_size, arma::vec t_accept) {
  
  // first, Uniform sampling on all categories of the parents
  arma::vec u(par_grid.n_rows);
  for (arma::uword x_par = 0; x_par < par_grid.n_rows; ++x_par) {
    u(x_par) = beta_update_rcpp(arma::sum(t), mct(x_par));
  }
  
  // next, MALA on all categories of the child
  for (int x = 0; x <= k; ++x) {
    // find indices where the first column of par_child_grid equals x
    arma::uvec indices = arma::find(par_child_grid.col(0) == x);
    
    // extract the corresponding N values from ct
    arma::vec N = ct.elem(indices);
    
    // perform MALA update
    Rcpp::List t_x_mala = mala_update_rcpp(t(x), step_size(x), b, rho, k, N, u);
    
    // update t and t_accept with the results
    t(x) = Rcpp::as<double>(t_x_mala["next_state"]);
    t_accept(x) = Rcpp::as<double>(t_x_mala["acceptance"]);
  }
  
  return Rcpp::List::create(Rcpp::Named("t_update") = t,
                            Rcpp::Named("u_update") = u,
                            Rcpp::Named("t_accept") = t_accept);
}


// [[Rcpp::export]]
Rcpp::List mala_edge_sampler_rcpp(int R,
                                  const arma::vec& x_j, const arma::mat& X_minus,
                                  double b_j, double rho_j,
                                  double c_j, double d_j,
                                  arma::vec t_init,
                                  const arma::vec& step_size, int stops = 1) {
  
  // Initialize
  int k_j = arma::max(x_j);
  int num_cols = X_minus.n_cols;
  
  // initialize z as an arma::uvec
  arma::uvec z(num_cols);
  double prob = c_j / (c_j + d_j);
  for (int i = 0; i < num_cols; ++i) {
    z(i) = static_cast<arma::uword>(R::rbinom(1, prob));
  }
  
  arma::umat z_samps(R, num_cols);
  z_samps.row(0) = z.t();
  
  arma::vec t = t_init;
  arma::mat t_samps(R, k_j + 1);
  t_samps.row(0) = t.t();
  
  arma::vec t_accept = arma::zeros(k_j + 1);
  arma::mat t_accept_samps(R, k_j + 1);
  t_accept_samps.row(0) = t_accept.t();
  
  // run Gibbs sampler
  for (int r = 1; r < R; ++r) {
    // update latent variables
    Rcpp::List current_parents;
    if (arma::sum(z) == 0) {
      current_parents = create_grids_and_table_root_rcpp(x_j);
    } else {
      // find indices where z is 1
      arma::uvec parent_indices = arma::find(z == 1);
      arma::mat x_Pa_j = X_minus.cols(parent_indices);
      current_parents = create_grids_and_table_child_rcpp(x_j, x_Pa_j);
    }
    
    Rcpp::List mwg_update = mala_within_gibbs_graph_rcpp(
      t,
      Rcpp::as<arma::vec>(current_parents["ct"]),
      Rcpp::as<arma::vec>(current_parents["mct"]),
      Rcpp::as<arma::mat>(current_parents["par_grid"]),
      Rcpp::as<arma::mat>(current_parents["par_child_grid"]),
      b_j, rho_j, k_j,
      step_size, arma::zeros(k_j + 1)
    );
    
    t = Rcpp::as<arma::vec>(mwg_update["t_update"]);
    t_accept = Rcpp::as<arma::vec>(mwg_update["t_accept"]);
    
    // update graph
    z = update_parents_rcpp(t, z, x_j, X_minus, c_j, d_j);
    
    // exporting samples
    t_samps.row(r) = t.t();
    t_accept_samps.row(r) = t_accept.t();
    z_samps.row(r) = z.t();
    
    // print stops
    if ((r + 1) % stops == 0) {
      Rcpp::Rcout << r + 1 << std::endl;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("t") = t_samps,
    Rcpp::Named("t_accept") = t_accept_samps,
    Rcpp::Named("z") = z_samps
  );
}


// [[Rcpp::export]]
double log_posterior(const arma::vec& t,
                     const arma::vec& z,
                     const arma::vec& x_j,
                     const arma::mat& X_minus,
                     double b_j,
                     double rho_j,
                     double c_j,
                     double d_j) {
  // prelims
  double beta_j = arma::sum(t);
  int k_j = x_j.max();
  int M = arma::sum(z);
  int W = arma::sum(1-z);
  
  // make table and grid
  Rcpp::List obj;
  int K_Pa_j;
  if (M == 0) {
    obj = create_grids_and_table_root_rcpp(x_j);
    K_Pa_j = 1;
  } else {
    obj = create_grids_and_table_child_rcpp(x_j, X_minus.cols(arma::find(z == 1)));
    K_Pa_j = arma::prod(arma::max(X_minus.cols(arma::find(z == 1)), 0)+1);
  }
  arma::vec ct = Rcpp::as<arma::vec>(obj["ct"]);
  arma::vec mct = Rcpp::as<arma::vec>(obj["mct"]);
  arma::mat par_grid = Rcpp::as<arma::mat>(obj["par_grid"]);
  arma::mat par_child_grid = Rcpp::as<arma::mat>(obj["par_child_grid"]);
  // first term
  double A = 0.0;
  for (arma::uword x = 0; x < par_grid.n_rows; ++x) {
    A += lgamma(beta_j) - lgamma(mct(x) + beta_j);
  }
  double B = 0.0;
  for (arma::uword x = 0; x <= k_j; ++x) {
    B += -b_j*t(x) + ((rho_j/k_j*1.0) - 1)*log(t(x))-K_Pa_j*lgamma(t(x));
  }
  double C = 0.0;
  for (arma::uword x = 0; x < par_child_grid.n_rows; ++x) {
    C += lgamma(ct(x) + t(par_child_grid(x, 0)));
  }
  double D = lgamma(M+c_j) + lgamma(W + d_j);
  return(A+B+C+D);
}

// [[Rcpp::export]]

arma::vec log_posterior_samples(const arma::mat& t,
                                const arma::mat& z,
                                const arma::vec& x_j,
                                const arma::mat& X_minus,
                                double b_j,
                                double rho_j,
                                double c_j,
                                double d_j) {
arma::vec lp_mcmc = arma::zeros(t.n_rows,1);
  for (int r=0; r < t.n_rows; ++r) {
    lp_mcmc(r) = log_posterior(t.row(r).t(),
                  z.row(r).t(),
                  x_j,
                  X_minus,
                  b_j,
                  rho_j,
                  c_j,
                  d_j);
  }
  return(lp_mcmc);
}

// [[Rcpp::export]]
int model_update_cpp(arma::vec t, List ct, List mct, List par_grid, List par_child_grid,
                     arma::vec prior_probs) {
  int M = ct.size();
  arma::vec log_model_probs = arma::zeros(M,1);
  for (int m = 0; m<M; m++) {
    log_model_probs[m] = log(prior_probs(m)) +
      log_node_score_rcpp(
      t,
      Rcpp::as<arma::vec>(ct[m]),
      Rcpp::as<arma::vec>(mct[m]),
      Rcpp::as<arma::mat>(par_grid[m]),
      Rcpp::as<arma::mat>(par_child_grid[m])
    );
  }
  arma::vec prbs = arma::exp(log_model_probs - arma::max(log_model_probs));
  return(sample_arma(prbs));
}

// [[Rcpp::export]]
arma::mat debug(List par_grid, int pset) {
  return Rcpp::as<arma::mat>(par_grid[pset - 1]);
}

// [[Rcpp::export]]
Rcpp::List mala_model_sampler_rcpp(int R,
                                   int k_j,
                                  const List& ct, const List& mct, 
                                  const List& par_grid, const List& par_child_grid,
                                  double b_j, double rho_j,
                                  arma::vec prior_probs,
                                  arma::vec t_init, int pset_init,
                                  const arma::vec& step_size, int stops = 1) {
  
  // Initialize
  int pset = pset_init;
  arma::vec pset_samps = arma::zeros(R,1);
  
  arma::vec t = t_init;
  arma::mat t_samps(R, k_j + 1);
  t_samps.row(0) = t.t();
  
  arma::vec t_accept = arma::zeros(k_j + 1);
  arma::mat t_accept_samps(R, k_j + 1);
  t_accept_samps.row(0) = t_accept.t();
  //Rcout << "Initiazliation Done";
  // Run Gibbs sampler
  for (int r = 1; r < R; ++r) {
    // update latent variables
    Rcpp::List mwg_update = mala_within_gibbs_graph_rcpp(
      t,
      Rcpp::as<arma::vec>(ct[pset-1]),
      Rcpp::as<arma::vec>(mct[pset-1]),
      Rcpp::as<arma::mat>(par_grid[pset-1]),
      Rcpp::as<arma::mat>(par_child_grid[pset-1]),
      b_j, rho_j, k_j,
      step_size, arma::zeros(k_j + 1)
    );
    //Rcout << "t sampled";
    
    t = Rcpp::as<arma::vec>(mwg_update["t_update"]);
    t_accept = Rcpp::as<arma::vec>(mwg_update["t_accept"]);
    
    // Update graph
    pset = model_update_cpp(t, ct, mct, par_grid, par_child_grid,
                            prior_probs);
    //Rcout << "Graph updated";
    
    // Exporting samples
    t_samps.row(r) = t.t();
    t_accept_samps.row(r) = t_accept.t();
    pset_samps(r) = pset;
    //Rcout << "Samples saved";
    
    // Print stops
    if ((r + 1) % stops == 0) {
      Rcpp::Rcout << r + 1 << std::endl;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("t") = t_samps,
    Rcpp::Named("t_accept") = t_accept_samps,
    Rcpp::Named("Pset") = pset_samps
  );
}
