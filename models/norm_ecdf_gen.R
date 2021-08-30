# Instead of reconstructing a new ecdf function for every run drawn from
# the same norm(mean = 0, sd = 1), use an standard ecdf averaged from multiple runs
# I suspect this would in the end remove the need to sample from a norm distribution
# cos it elminates the norm variation, but I'd rather run this once to confirm
# -Fai

# Average ecdf based on method written here
# https://stats.stackexchange.com/questions/4261/should-i-use-an-average-ecdf
# Our issue of local maxima and minima discussed in link should not matter because 
# we have have a normalized clear range of 0-1

# Values to impute
repeats <- 1e3
impute <- seq(0, 1, length = 1e3)
ecdfs <- matrix(NA, ncol = repeats, nrow = length(impute))

for(i in 1:repeats){
    # Copied from original code
    probs <- rnorm(1e5) #randomly draw 100 probabilities from a normal distribution
    probs.normalized <- (probs - min(probs))/(max(probs) - min(probs)) #shift distribution to start at 0 and have 
    p <- ecdf(probs.normalized) #cumulative distribution
    ecdfs[,i] <- p(impute)
    print(i)
}

mean_ecdf = rowMeans(ecdfs)
norm_ecdf <- stepfun(impute, c(0, mean_ecdf))
plot(norm_ecdf)

save(norm_ecdf,file = "norm_ecdf.Rdata")



