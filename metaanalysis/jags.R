
###############################################
## Model
###############################################

jags.mod <- function(){
  
  # Likelihood:
  for (i in binom_dist){
    
    n_ind_outcome[i] ~ dbin(mu[i], n_ind_contributedsamples[i])
    
    logit(mu[i]) <- a[study_id[i]] + # intercept - each trial has an intercept
      b[study_id[i]] * abx_dur[i]  + # each trial has a slope
      c[study_id[i]] * time_baseline_endoffu[i]
    
    #loglike[i] <- dbin(outcome_incub5[i], mu[i], 1) # For WAIC computation
    
  }
  
  for (i in pois_dist){
    
    n_ind_outcome[i] ~ dpois(mu[i] * n_ind_contributedsamples[i])
    
    logit(mu[i]) <- a[study_id[i]] +   # intercept - each trial has an intercept
      b[study_id[i]] * abx_dur[i]      # each trial has a slope
  }
  
  for (w in 1:study_N){
    
    a[w] ~ dnorm(a0, sigma_a);
    b[w] ~ dnorm(b0, sigma_b);
    c[w] ~ dnorm(c0, sigma_c);
    
    # a[w] <- a0 + aprimed[w] * sigma.a
    # b[w] <- b0 + bprimed[w] * sigma.b
    # c[w] <- c0 + cprimed[w] * sigma.c
    # 
    # aprimed[w] ~ dnorm(0, 5);T(0,);
    # bprimed[w] ~ dnorm(0, 5);T(0,);
    # cprimed[w] ~ dnorm(0, 5);T(0,);
    
  }
  
  # Priors:
  # a0 ~ dnorm(0.1, 1);T(0,);
  # sigma.a ~ dunif(0, 5);
  # b0 ~ dnorm(0.1, 1);
  # sigma.b ~ dunif(0, 5);
  # c0 ~ dnorm(0.1, 1);
  # sigma.c ~ dunif(0, 5);
  
  a0 ~ dnorm(1, 5);T(0,);
  b0 ~ dnorm(0, 5);
  c0 ~ dnorm(0, 5);T(0,);
  
  sigma_a ~ dunif(0, 2); 
  sigma_b ~ dunif(0, 2); 
  sigma_c ~ dunif(0, 2);
  
  
}

init <- function(){
  list(a0 = 0.1, b0 = 0.1, c0 = 0.1, 
       sigma_a = 1, sigma_b = 1, sigma_c = 1)
}

params <- c("a", "b", "c", "a0", "b0", "sigma_a", "sigma_b")

