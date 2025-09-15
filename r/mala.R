log_h <- function(t, b, rho, k, N, U) {
  # un-normalized log of h
  g <- rho/(k+1)
  return( -b*t + (g-1)*log(t) + sum(lgamma(N+t) - lgamma(t) + t*log(U)))
}

grad_log_h <- function(t, b, rho, k, N, U) {
  # gradient of the logarithm of h
  g <- rho/(k+1)
  return( - b + ((g - 1)/t) + sum(digamma(N+t) - digamma(t) + log(U)) )
}

mala_proposal <- function(t_current, step_size, 
                              b, rho, k, N, U) {
  # langevin proposal for MALA
  # assumes that the mass matrix is the identity
  gradient <- grad_log_h(t_current, b, rho, k, N, U)
  t_proposal <- t_current + (step_size^2/2) * gradient +
    step_size * rnorm(1) # Assuming identity mass matrix
  
  return(t_proposal)
}

log_mala_density <- function(t_current, t_proposal, step_size,
                             b, rho, k, N, U) {
  # log Gaussian density for langevin proposal, needed for acceptance prob
  return(dnorm(x=t_proposal,
               mean=t_current + (step_size^2/2)*grad_log_h(t_current, b, rho, k, N, U),
               sd=step_size,
               log=TRUE))
}

mala_update <- function(t_current, step_size,
                        b, rho, k, N, U) {
  # compute proposal
  t_proposal <- mala_proposal(t_current, step_size, 
                                  b, rho, k, N, U)
  if (t_proposal<0) {
    return(list(next_state = t_current, acceptance = 0))
  } else {
    log_acceptance_ratio <- log_h(t_proposal, b, rho, k, N, U) - log_h(t_current, b, rho, k, N, U) +
      log_mala_density(t_proposal=t_current, t_current=t_proposal, step_size, b, rho, k, N, U) -
      log_mala_density(t_proposal=t_proposal, t_current=t_current, step_size, b, rho, k, N, U)
    # acceptance probability
    if (log(runif(1)) < log_acceptance_ratio) {
      return(list(next_state = t_proposal, acceptance = 1))
    } else {
      return(list(next_state = t_current, acceptance = 0))
    }
  }
}

### sampling the Us's ###
beta_update <- function(beta, n) {
  return(rbeta(1,beta,n))
}

