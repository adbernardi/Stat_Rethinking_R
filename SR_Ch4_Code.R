# Ch 4 Code - Linear Models

pos <- replicate(1000, sum(runif(16,-1,1)))
plot(density(pos))
hist(pos)

?runif

# further investigation normal distribution with the loci code 
prod(1 + runif(12,0,0.1))
?prod

growth <- replicate(10000, prod(1 + runif(12,0,0.1)))
library(rethinking)
dens(growth, norm.comp=TRUE)

# checking how big and small numbers work with the normal distribution
big <- replicate(10000, prod(1 + runif(12,0,0.5)))
small <- replicate(10000, prod(1 + runif(12,0,0.1)))

dens(small, norm.comp = TRUE)
dens(big, norm.comp = TRUE)

# smaller numbers converge, multiplying is the same as adding for small #'s

# log scaling for larger numbers 

log.big <- replicate(10000, log(prod(1 + runif(12,0,0.5))))
dens(log.big, norm.comp = TRUE)

dnorm(0,0,0.1)
# the above is a rate change of the probability, the probability density 

# 4.2 using the lexicon and language of modeling 
# from earlier, using the binomial distribution and uniform prior in the globe tossing example 

# calcluating in a way that'll look familiar 
w <- 6 # number of observed water tosses 
n <- 9 # number of trials 
p_grid <- seq(from=0, to=1, length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)
posterior

# looks familiar 
# Eventually building a Gaussian regression model using Bayesian updating 

# 4.3.1 - starting with the data 
# loading
library(rethinking)
data(Howell1)
d <- Howell1

# snooping around 
str(d)

# working with just the height for now 
d$height

# dealing only with adults, for simplicity of height 
head(d)
d2 <- d[d$age >= 18,]
d2
str(d2)
d2$height
head(d2)
# doing some visualizing 
dens(d2$height)

# starting with picking priors of the parameters we're going to estimate 
# starting with the prior for mean height and deviation
curve(dnorm(x, 178, 20), from=100, to=250)
curve(dnorm(x, 178, 0.1), from=100, to=250)
# and a flat prior for SD 
curve(dunif(x, 0, 50), from=10, to=60)
# density of a uniform dist.

# drawing samples from these distributions
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)

prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# sandboxing
# a more aggressive example w a smaller spread?
#sample_mu <- rnorm(1e4, 178, 15)
#prior_h <- rnorm(1e4, sample_mu, sample_sigma)
#dens(prior_h)

#sample_mu <- rnorm(1e4, 178, 30)
#prior_h <- rnorm(1e4, sample_mu, sample_sigma)
#dens(prior_h)

# going back to the original 
# 4.3.3. Grid approximation of the posterior distribution
# putting in the code for this computation, returning to how it works specifically later

mu.list <- seq(from=140, to=160, length.out=200 )
sigma.list <- seq(from=4, to=9, length.out=200)
post <- expand.grid(mu=mu.list, sigma=sigma.list)
post$LL <- sapply(1:nrow(post), function(i) sum(dnorm(d2$height, 
                                                      mean=post$mu[i],
                                                      sd=post$sigma[i],
                                                      log=TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) + dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))

#contour plot
contour_xyz(post$mu, post$sigma, post$prob)
# heatmap 
image_xyz(post$mu, post$sigma, post$prob)

# sampling from the posterior distribution as before, but now with two parameters 
sample.rows <- sample(1:nrow(post), size=1e4, replace=TRUE , prob=post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

# now plotting
plot(sample.mu, sample.sigma, cex=0.5, pch=16, col=col.alpha(rangi2, 0.1))

# now to think and describe the samples 
dens(sample.mu)
dens(sample.sigma)

# summarizing the widths of these densities using a familiar tool by now 
HPDI(sample.mu)
HPDI(sample.sigma)

# sampling 20 randm heights to illustrate the potential pitfall with sigma 
d3 <- sample(d2$height, size=20)

# this is going to be the same as the code from above 
mu.list <- seq(from=150, to=170, length.out=200)
sigma.list <- seq(from=4, to=20, length.out=200)
post2 <- expand.grid(mu=mu.list, sigma=sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i) 
    sum(dnorm(d3, mean=post2$mu[i], sd=post2$sigma[i], log=TRUE)))
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, TRUE) + (dunif(post2$sigma, 0, 50, TRUE))
post2$prob <- exp(post2$prod - max(post2$prod))
sample2.rows <- sample(1:nrow(post2), size=1e4, replace=TRUE, prob=post2$prob)
sample2.mu <- post2$mu[sample2.rows]
sample2.sigma <- post2$sigma[sample2.rows]
plot(sample2.mu, sample2.sigma, cex=0.5, 
     col=col.alpha(rangi2, 0.1), 
     xlab="mu", ylab="sigma", pch=16)

# visualizing another way 
dens(sample2.sigma, norm.comp=TRUE)
# here we see the long positive tail

# 4.3.5. Fitting a model with map 
# making sure the data is loaded 

str(d)
str(d2)
fivenum(d2$age)
fivenum(d$age)

# approximating and translating the language of the model 
flist <- alist(
          height ~ dnorm(mu, sigma) ,
          mu ~ dnorm(178, 20) , 
          sigma ~ dunif(0,50)
)
flist
?alist()
# commas separate each line of the model definition 

?map()
m4.1 <- map(flist, data=d2)

precis(m4.1)

?precis()
# provides gaussian approximations for each parameter's marginal dist. 
# trying again after rebuilding the prior and its std. dev.
flist <- alist(
  height ~ dnorm(mu, sigma) ,
  mu ~ dnorm(178, 0.1) , 
  sigma ~ dunif(0,50)
)

m4.2 <- map(flist, data=d2)

precis(m4.2)

vcov(m4.1)
?vcov

# decomposing the variance/covariance matrix 

diag(vcov(m4.1))
cov2cor(vcov(m4.1))

# now sampling from the multi-dimensional gaussian dist 
post <- extract.samples(m4.1, n=1e4)
head(post)

precis(post)

?precis

m4.1_logsigma <- map(
                  alist(
                      height ~ dnorm(mu, exp(log_sigma)) , 
                      mu ~ dnorm(178, 20) , 
                      log_sigma ~ dnorm(2, 10)
                  ) , data=d2)


# getting back on the natural scale 
post <- extract.samples(m4.1_logsigma)
sigma <- exp(post$log_sigma)

# 4.4 Adding a predictor

# adding a predictor to get the whole model feel 
# starting with some visual EDA 

plot(d2$height ~ d2$weight)

# 4.4.2 fitting the model, now with a predictor
library(rethi)

# loading the data again 
# now fitting the model 
m4.3 <- map(
      alist(
          height ~ dnorm(mu, sigma) ,
          mu <- a + b*weight ,
          a ~ dnorm(178, 100) ,
          b ~ dnorm(0, 10) ,
          sigma ~ dunif(0, 50)
      ),
    data =d2)

# inspecting the estimates 
precis(m4.3)

# let's now move on to the correlation matrix 
# seeing how our parameters relate to one another 
cov2cor(vcov(m4.3))

# another way of doing this 
precis(m4.3, corr = TRUE)

# tricks to avoid the strong negative correlation b/w our alpha and beta 
# starting with centering, subtracting the mean of a variable from each value 

d2$weight.c <- d2$weight - mean(d2$weight)

mean(d2$weight.c)

# refitting the model with this center and seeing what it gives us 
m4.4 <- map(
  alist(
    height ~ dnorm(mu, sigma) ,
    mu <- a + b*weight.c ,
    a ~ dnorm(178, 100) ,
    b ~ dnorm(0, 10) ,
    sigma ~ dunif(0, 50)
  ),
  data =d2)

precis(m4.4)
cov2cor(vcov(m4.4))

precis(m4.4, corr=TRUE)

?map()

# now, we are going to superimpose the MAP values over existing plots of the height and weight data 

plot(height ~ weight, data=d2, pch=21, bg="blue")
?plot()
abline(a=coef(m4.3)["a"], b=coef(m4.3)["b"])

# adding uncertainty into our map line 

post <- extract.samples(m4.3)

post[1:5,]

# the following code pulls 10 samples and incorporates that into our plotting

N <- 10
dN <- d2[1:N, ]
mN <- map(
  alist(
    height ~ dnorm(mu, sigma) ,
    mu <- a + b*weight , 
    a ~ dnorm(178, 100) , 
    b ~ dnorm(0, 10) , 
    sigma ~ dunif(0, 50)
  ), data=dN)

# extracting and plotting 20 of these lines to get an idea 

post <- extract.samples(mN, n=20)

# displaying these data 

plot( dN$weight, dN$height, 
      xlim=range(d2$weight), ylim=range(d2$height), 
      col=rangi2, xlab="weight", ylab="height")
mtext(concat("N = ", N))

# now adding the lines 

for (i in 1:20 )
  abline(a=post$a[i] , b=post$b[i] , col=col.alpha("black", 0.3))

# trying to get the cloud around the regression line 
# this occurs by sampling both of the parameters at distinct times 

# taking an example of 50 for weight and sampling while holding 50 constant 

mu_at_50 <- post$a + post$b*50
head(mu_at_50)

dens(mu_at_50)

# trying another way of plotting and conveying this information 

dens(mu_at_50, col=rangi2, lwd=2, xlab="mu|weight=50")

HPDI(mu_at_50, prob=0.89)

# now we need to do this for all parameter combinations, not just 50

mu <- link(m4.3)
str(mu)

# now we want to do this for ALL weight values 
# define sequence of weights to compute predictions for 
# these values will be on the horizontal axis 
weight.seq <- seq(from=25, to=70, by=1)

# using link to compute mu for each sample from the posterior and for each weight in the above sequence (weight.seq)

mu <- link(m4.3, data=data.frame(weight=weight.seq))
str(mu)

# and now using type="n" to hide the raw data 
plot(height ~ weight, d2, type="n")

# looping over samples and plotting each mu value 
for (i in 1:1000) 
  points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2, 0.1))

# summarizing the distribution of mu 
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

str(mu.mean)
str(mu.HPDI)

# plotting raw data 
# fading points for visibility 
plot(height ~ weight, data=d2, col=col.alpha(rangi2, 0.5))

# adding the MAP line (mean mu for each weight)
lines(weight.seq, mu.mean)

# adding a shaded region to denote 89% HDPI
shade(mu.HPDI, weight.seq)

# now are going to move into predicting actual heights 
sim.height <- sim(m4.3, data=list(weight=weight.seq))

str(sim.height)

# summarizing 
height.PI <- apply(sim.height, 2, PI, prob=0.89)

# now plotting everything we have 

# raw data 
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))

# MAP line 
lines(weight.seq, mu.mean)

# HPDI region for line 
shade(mu.HPDI, weight.seq)

# draw PI region for simulated heights 
shade(height.PI, weight.seq)

# 4.5 Polynomial Regression 
# nothing special about straight lines!

data(Howell1)
d <- Howell1
str(d)

plot(height ~ weight, data=d)

# standardizing the weight variable 
d$weight.s <- (d$weight - mean(d$weight)) / sd(d$weight)

# checking 
plot(height ~ weight.s, data=d)

# now building a quadratic model 

d$weight.s2 <- d$weight.s^2
m4.5 <- map(
  alist(
    height ~ dnorm(mu, sigma) , 
    mu <- a + b1*weight.s + b2*weight.s2 , 
    a ~ dnorm(178, 100) , 
    b1 ~ dnorm(0, 10) , 
    b2 ~ dnorm(0, 10) , 
    sigma ~ dunif(0, 50)
  ) ,
data = d)

precis(m4.5)

# now plotting what we've found 

weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight.s=weight.seq, weight.s2=weight.seq^2)
mu <- link(m4.5, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

plot(height ~ weight.s, d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)


# Chapter 4 Exercises 

# Medium
# simulating observed heights from the prior 

?sim()

# translating the model to start with 
library(rethinking)

n <- 1000

mu <- rnorm(n, 0, 10)
sigma <- runif(n, 0, 10)

y <- rnorm(n, mu, sigma)

dens(y)
# 4M2 - Translating the model to a map formula 
# y ~ normal(mu, sigma)
# mu ~ normal(0, 10)
# sigma ~ uniform(0, 10)

# need ot actually load in the data?
?dnorm
?rnorm

# 4M3 - Translating the map model into a mathematical definition 

# done off screen 

# M4, M5, M6 all done off screen 

# 4H1 

# filling in a table with model-based predictions 

data("Howell1")
head(Howell1)

d <- Howell1

str(d)

# now time to set up a model framework

plot(height ~ weight, data = d)

# linear model probably won't cut it, quadratic it is
# we want to use weight to predict expected height 

# translating the model framework and trying to use this for prediction
d$weight.sq <- (d$weight)^2

h1_pred <- map(
  alist(
    height ~ dnorm(mu, sigma) , 
    mu ~ a + b1*weight + b2*weight.sq , 
    a ~ dnorm(150, 50) , 
    b1 ~ dnorm(0, 10) , 
    b2 ~ dnorm(0, 10) , 
    sigma ~ dunif(0, 50)
  ) , 
data = d)

# let's keep going with trying this 
precis(h1_pred)

# maybe weight does have to be standardized 

d$weight_std <- (d$weight - mean(d$weight)) / sd(d$weight)
d$weight_std_sq <- (d$weight_std)^2

h1_pred <- map(
  alist(
    height ~ dnorm(mu, sigma) , 
    mu ~ a + b1*weight_std + b2*weight_std_sq , 
    a ~ dnorm(150, 50) , 
    b1 ~ dnorm(0, 10) , 
    b2 ~ dnorm(0, 10) , 
    sigma ~ dunif(0, 50)
  ) , 
  data = d)

precis(h1_pred)

# seems better 
# we could even try with a cubic and see how that does 

# now trying to plot this 

weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight_std=weight.seq, weight_std_sq = (weight.seq)^2)
mu <- link(h1_pred, data=pred_dat)

pred_dat

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(h1_pred, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

# and now on to the actual plotting

plot(height ~ weight_std, d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

# a cubic would probably be more accurate, given people that are huge won't all of the sudden get shorter?

mean(d$weight)

# try a test case 
(46.95 - mean(d$weight)) / sd(d$weight)

precis(h1_pred)

# now to try and make an actual prediction
# starting with our array of weights 

weights <- c(46.95, 43.72, 64.78, 32.59, 54.63)

