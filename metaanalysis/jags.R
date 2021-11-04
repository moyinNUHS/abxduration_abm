
###############################################
## Model 1 with random effect 
###############################################

jags.rand <- function(){
  
  # Likelihood:
  for (i in 1:length(study_id)){ # for each arm 
    
    n_ind_outcome[i] ~ dbin(mu[i], n_ind_contributedsamples[i])
    
    logit(mu[i]) <- a[study_id[i]] + # intercept - each trial has an intercept
      b[study_id[i]] * abx_dur[i] + 
      c[study_id[i]] * time_baseline_endoffu[i] + 
      d[study_id[i]] * setting[i] 
    
    loglike[i] <- dbin(n_ind_outcome[i], mu[i], n_ind_contributedsamples[i]) # For WAIC computation
  }
  
  for (w in 1:study_N){ # for each study 
    
    a[w] ~ dnorm(a0, sigma_a);
    b[w] ~ dnorm(b0, sigma_b);
    c[w] ~ dnorm(c0, sigma_c);
    d[w] ~ dnorm(d0, sigma_d);
    
  }
  
  a0 ~ dnorm(mean_prior, mean_sd_prior);
  b0 ~ dnorm(mean_prior, mean_sd_prior); 
  c0 ~ dnorm(mean_prior, mean_sd_prior); 
  d0 ~ dnorm(mean_prior, mean_sd_prior); 
  
  sigma_a ~ dgamma(sigma_prior1, sigma_prior2); 
  sigma_b ~ dgamma(sigma_prior1, sigma_prior2); 
  sigma_c ~ dgamma(sigma_prior1, sigma_prior2); 
  sigma_d ~ dgamma(sigma_prior1, sigma_prior2);
  
}

init.rand <- function(){
  list(a0 = 0, b0 = 0,  c0 = 0, d0 = 0, sigma_a = 1, sigma_b = 1, sigma_c = 1, sigma_d = 1)
}

params.rand <- c("a", "b", "c", "d", 
                 "a0", "b0", "c0", "d0", 
                 "sigma_a", "sigma_b",  "sigma_c", "sigma_d",
                 "loglike")

###############################################
## Model 2 with random effect, no D
###############################################

jags.randnod <- function(){
  
  # Likelihood:
  for (i in 1:length(study_id)){ # for each arm 
    
    n_ind_outcome[i] ~ dbin(mu[i], n_ind_contributedsamples[i])
    
    logit(mu[i]) <- a[study_id[i]] + # intercept - each trial has an intercept
      b[study_id[i]] * abx_dur[i] + 
      c[study_id[i]] * time_baseline_endoffu[i] 
    
    loglike[i] <- dbin(n_ind_outcome[i], mu[i], n_ind_contributedsamples[i]) # For WAIC computation
  }
  
  for (w in 1:study_N){ # for each study 
    
    a[w] ~ dnorm(a0, sigma_a);
    b[w] ~ dnorm(b0, sigma_b);
    c[w] ~ dnorm(c0, sigma_c);
    
  }
  
  a0 ~ dnorm(mean_prior, mean_sd_prior);
  b0 ~ dnorm(mean_prior, mean_sd_prior); 
  c0 ~ dnorm(mean_prior, mean_sd_prior); 
  
  sigma_a ~ dgamma(sigma_prior1, sigma_prior2); 
  sigma_b ~ dgamma(sigma_prior1, sigma_prior2); 
  sigma_c ~ dgamma(sigma_prior1, sigma_prior2); 
  
}

init.randnod <- function(){
  list(a0 = 0, b0 = 0,  c0 = 0, sigma_a = 1, sigma_b = 1, sigma_c = 1)
}

params.randnod <- c("a", "b", "c",
                 "a0", "b0", "c0",
                 "sigma_a", "sigma_b",  "sigma_c",
                 "loglike")


###############################################
## Model 3 with no random effect 
###############################################

jags.norand <- function(){
  
  # Likelihood:
  for (i in 1:length(study_id)){ # for each arm 
    
    n_ind_outcome[i] ~ dbin(mu[i], n_ind_contributedsamples[i])
    
    logit(mu[i]) <- a + # intercept - each trial has an intercept
      b * abx_dur[i] + 
      c * time_baseline_endoffu[i] + 
      d * setting[i] 
      
      loglike[i] <- dbin(n_ind_outcome[i], mu[i], n_ind_contributedsamples[i]) # For WAIC computation
  }
  
  a ~ dnorm(mean_prior, mean_sd_prior);
  b ~ dnorm(mean_prior, mean_sd_prior);
  c ~ dnorm(mean_prior, mean_sd_prior);
  d ~ dnorm(mean_prior, mean_sd_prior);
  
}


init.norand <- function(){
  list(a = 0, b = 0,  c = 0, d = 0)
}

params.norand <- c("a", "b", "c", "d", "loglike")

