# written by K. Garner 20th Sept, 2017
# code runs through equations # of Wagenmakers, EJ (2007) A practical solution to
# the pervasive problems of p-values. Psychonomic Bulletin and Review

# if you don't have the wesanderon package installed, please install it now for the figures to work
# install.packages("wesanderson")
rm(list = ls())

#### SECTION: Bayesian parameter estimation
# ---------------------------------------------------------------------------------------
# eq (1) 
# binomial probability of getting 9/12 answers correct on true/false questions by chance
pDgtheta = dbinom(9, size = 12, prob = 0.5)

# now to achieve equation 5 (Pr(theta|D))
# estimating the prior
# generate a uniform (uninformed prior)
x = seq(0, 1, 0.001) # possible theta values
un_prior = dbeta(x,1,1) # uninformative prior - using the beta function to generate density
par(mfrow=c(3,1), mar=c(2,5,2,5))
plot(x, un_prior, main = "Beta(alpha = 1, beta = 1) Prior Distribution", type = "l", lwd = 2,
     ylab = "density", xlab = "binomial parameter - theta",
     col = wesanderson::wes_palette("FantasticFox")[3])

# PrD = 1/(n+1) - we have an analytic solution to PrD because with this simple form, we can use
# the gamma function to approximate the binomial model for PrD

# so Pr(theta|D) = (n+1)*Pr(s|theta,n)
# where s = 9, n = 12, theta = x
Pr_theta_g_D = 13*choose(12,9)*x^9*(1-x)^3
plot(x, Pr_theta_g_D, main = "Pr(theta|d) - aka posterior distribution (binom model)", type = "l", lwd = 2,
     ylab = "density", xlab = "binomial parameter - theta",     
     col = wesanderson::wes_palette("FantasticFox")[1])

# when updating a Beta prior through a binomial likelihood function, the resulting posterior 
# distribution is still a beta distribution (Beta(alpha + s, beta + n - s)) # sorry, lower case beta
# refers to the greek letter - which may be somewhat confusing
beta_posterior = dbeta(x,10,4) #(9 + 1), (3 + 1) - the 1s refer to the prior
plot(x, beta_posterior, main = "Pr(theta|d) - aka posterior distribution (Beta(10,4))", type = "l", lwd = 2,
     ylab = "density", xlab = "binomial parameter - theta",     
     col = wesanderson::wes_palette("FantasticFox")[5])

#### SECTION: Bayesian Hypothesis Testing
# ---------------------------------------------------------------------------------------
# visualising priors H0 and H1
# H0 = .5 (chance)
par(mfrow = c(2,1), mar = c(2,5,2,5))
Pr_H0 = rep(0, length(x))
Pr_H0[x == 0.5] = 1
plot(x, Pr_H0, main = "H0: theta = .5", type = "l", lwd = 2, 
     ylab = "density", xlab = "binomial parameter - theta",
     col = wesanderson::wes_palette("FantasticFox")[3])

# Pr H1 (this is the same as above)
plot(x, un_prior, main = "H1: Beta(alpha = 1, beta = 1)", type = "l", lwd = 2,
     ylab = "density", xlab = "binomial parameter - theta",
     col = wesanderson::wes_palette("FantasticFox")[3])
Pr_theta_g_H0 = choose(12,9)*0.5^12

#### SECTION: A Bayesian test of the p Postulate
# ---------------------------------------------------------------------------------------
n = seq(50, 10000, 50) # n's 50 to 10000
# for each level, find the # of obs that yields p ~ .05 using a binomial test (p(data|theta))
# H0 = 0.5
H0 = 0.5
get.obs.for.p <- function(n, H0){
  # this function works out how many observations required to get ~ p = .05 given the null hyp (theta = .5) is true
  # outcomes match those reported by Wagenmakers 2007 (at least the outcomes reported)
  upper_obs = seq(n/2,n,1)
  upper_ps = pbinom(upper_obs, n, prob = H0)
  idx = which(upper_ps-.975 > 0)[1]
  x = list(p = 1-upper_ps[idx], ob = upper_obs[idx]) 
}
NHST_results = lapply(n, FUN = get.obs.for.p, H0) # returns the # of observations per sample, and the p value 

# generate 2 Beta priors
x = seq(0, 1, 0.01) # possible theta values
un_prior = dbeta(x,1,1) # uniform prior
inf_prior = dbeta(x,6,3) # informed prior

# now for each observation/n, compute the posterior probabilities
# compute bayes factor (= ratio of predictive probabilities)
mat = matrix(unlist(NHST_results),ncol=2,byrow=TRUE)
s = mat[,2] # samples (correct beer guesses)
n = n # just to remind me (total beer guesses)
ps = mat[,1]*2
# uniform prior 
get.pH0.unf.prior <- function(s, n){
  # for each n - taking equation 8 from the paper
  H0 = .5
  BF = (n+1)*(dbinom(s, n,  prob = H0))
  PrH0_D = BF/(1+BF)
  return(PrH0_D)
}
pH0.unf.prior = mapply(get.pH0.unf.prior, s, n)

# inf prior
get.pH0.inf.prior <- function(s, n){
  # H1 stuff
  x = seq(0, 1, 0.01)
  inf_alpha = 6 # values for informed priors
  inf_beta = 3 # values for informed priors
  inf_prior = dbeta(x, inf_alpha, inf_beta)
  # Pr(D|H1) = 
  pD_H1 = sum((dbinom(s, n, prob = x)*inf_prior)/length(x)) # this is a rough numerical integration
  # then pH0
  H0 = 0.5
  BF = dbinom(s, n,  prob = H0)/pD_H1
  PrH0_D = BF/(1+BF)
  
  return(PrH0_D)
}
pH0.inf.prior = mapply(get.pH0.inf.prior, s, n)

# now use maximum likely theta as the prior
get.pH0.oracle.point.prior <- function(s, n){
  x = seq(0, 1, 0.01)
  pD_H1 = max(dbinom(s, n, x))
  BF = dbinom(s, n,  prob = H0)/pD_H1
  PrH0_D = BF/(1+BF)
  
  return(PrH0_D)
}
pH0.oracle.point.prior = mapply(get.pH0.oracle.point.prior, s, n)

# plot it all!
par(mfrow=c(2,1), mar=c(2,5,2,5))
# priors
plot(x, un_prior, main = "Two Priors", type = "l", lwd = 4,
     ylab = "density", ylim = c(0,3), xlab = "binomial parameter",
     col = wesanderson::wes_palette("FantasticFox")[1])
points(x, inf_prior, type = "l", lwd = 4,
       col = wesanderson::wes_palette("FantasticFox")[3])

# posterior odds
plot(n, pH0.unf.prior, type = "l", lwd = 3,
     ylab = "Posterior probability for H0", ylim = c(0,1), xlab = "n",
     col = wesanderson::wes_palette("FantasticFox")[1])
points(n, pH0.inf.prior, type = "l", lwd = 3,
       col = wesanderson::wes_palette("FantasticFox")[2])
points(n, pH0.oracle.point.prior, type = "l", lty = 2, lwd = 3,
       col = wesanderson::wes_palette("FantasticFox")[3])
par(new = T)
plot(n, ps, axes = F, ylim = c(0,1), ylab = "", xlab = "n",
     type = "l", lwd = 3,
     col = wesanderson::wes_palette("FantasticFox")[4])
axis(4, ylim=c(0,1), lwd=2)
mtext(4,text="NHST p value", line = 3)
legend(4000, 0.6, legend = c("Beta(1,1) prior", "Beta(6,3) prior", "Oracle", "NHST p-value"), lty = c(1,1,2,1),
       lwd = 3,
       col = c(wesanderson::wes_palette("FantasticFox")[1], wesanderson::wes_palette("FantasticFox")[2],
               wesanderson::wes_palette("FantasticFox")[3], wesanderson::wes_palette("FantasticFox")[4]))






