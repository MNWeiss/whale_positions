require(runjags) # you'll need to install these before running
require(rjags)

N <- 30 # number of individuals
G <- 100 # number of groups
p_present <- 0.5 # probability each individual is present in each group

# effects of position
b1 <- -5
b2 <- 5

par(mar = c(4,4,1,1))
layout(matrix(c(1,1,2,3),nrow=2))

curve(exp(b1*x + b2*(x^2)), xlab = "Lateral Position", ylab = "Relative Leadership Likelihood")

# individual effect SD
sigma_ind <- 0.5

curve(dnorm(x, sd = sigma_ind), xlim = c(-2,2), xlab = "Individual Effect", ylab = "Density")

# AR effects
lambda <- 2 # the decay coefficient for the autocorrelation
ar_beta_0 <- 2 # the baseline effect of autocorrelation

curve(ar_beta_0*exp(-lambda*x), xlab = "Time Since Last Observation", ylab = "Effect of Prior Leadership")

par(mfrow = c(1,1))

# Simulate data

# Data matrices
# In both matrices, rows are groups (or observations), and columns are individuals
# In present, a value 1 means the individual was present in that group, and 0 means they were absent
# In was_leader, a value of 1 means the individual was the leader in the most recent group, and 0 means they weren't
# In position, individuals in the group get a value from 0 to 1 indicating their lateral position (left to right, say)

present <- was_leader <- position <- matrix(0, nrow = G, ncol = N)

# Leadership data is just a vector, giving the index of the leader in each group

leader <- NA

# Lag is the time lag between each group; the scale here is from 0 to 1 but it doesn't matter in practice

lag <- runif(G)

# Individual random effects, which could be thought of as leadership tendency

ind_RE <- rnorm(n = N, mean = 0, sd = sigma_ind)

for(g in 1:G){
  
  present[g,] <- rbinom(n = N, size = 1, prob = p_present) # determine which individuals are present
  position[g,present[g,]==1] <- (sample(sum(present[g,]))-1)/(sum(present[g,])-1) # assign positions at random
  if(g == 1){
    was_leader[g,] <- 0 # in first group, no previous leadership data
  }else{
    was_leader[g,] <- ifelse(1:N == leader[(g-1)],1,0) # otherwise, record previous leader
  }
  
  ar_beta <- exp(-lambda*lag[g])*ar_beta_0 # determine effect of previous leadership given time lag
  
  p_leader <- exp(b1*position[g,] + b2*(position[g,]^2) + ind_RE + was_leader[g,]*ar_beta) # get leadership weights
  leader[g] <- sample(which(present[g,] == 1), prob = p_leader[present[g,]==1], 1) # sample leader from present whales, with weights
  
}
  
# Stick data together in a list so JAGS can use it

leadership_data <- list(
  N = N,
  G = G,
  leader = leader,
  present = present,
  position = position,
  lag = lag,
  was_leader = was_leader
)

# estimate the model from simulated data
leadership_fit <- run.jags("position model.jags",
                           data = leadership_data,
                           monitor = c("b1","b2","lambda","ar_beta_0","sigma_ind"))
