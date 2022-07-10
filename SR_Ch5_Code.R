# Chapter 5 - Multivariate Linear Models 

# loading in our data 

library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

# standardizing the predictor variable

d$MedianAgeMarriage.s <- (d$MedianAgeMarriage-mean(d$MedianAgeMarriage)) / sd(d$MedianAgeMarriage)

# fit the model 

m5.1 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma) , 
    mu <- a + bA * MedianAgeMarriage.s , 
    a ~ dnorm(10, 10) , 
    bA ~ dnorm(0, 1) , 
    sigma ~ dunif(0, 10)
  ), data = d
)

# computing the percentile interval of the mean 

MAM.seq <- seq(from=-3, to = 3.5, length.out = 30)
mu <- link(m5.1, data=data.frame(MedianAgeMarriage.s = MAM.seq))
mu.PI <- apply(mu, 2, PI)

# plotting all this 

plot(Divorce ~ MedianAgeMarriage.s, data = d, col=rangi2)
abline(m5.1)
shade(mu.PI, MAM.seq)

?link

# doing this for another plot 

d$Marriage.s <- (d$Marriage - mean(d$Marriage)) / sd(d$Marriage)
m5.2 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma) , 
    mu <- a + bR * Marriage.s , 
    a ~ dnorm(10, 10) , 
    bR ~ dnorm(0, 1) , 
    sigma ~ dunif(0, 10)
  ), data = d
)

# now plotting this 

plot(Divorce ~ Marriage.s, data = d, col=rangi2)
abline(m5.2)
shade(mu.PI, MAM.seq)

# trying to nail the shading 

#MAM.seq <- seq(from=-3, to = 3.5, length.out = 30)
mu <- link(m5.2, data=data.frame(Marriage.s = MAM.seq))
mu.PI <- apply(mu, 2, PI)


# now we'll use the multivariable model 

m5.3 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma) , 
    mu <- a + bR*Marriage.s + bA*MedianAgeMarriage.s , 
    a ~ dnorm(10, 10) , 
    bR ~ dnorm(0, 1) , 
    bA ~ dnorm(0, 1) ,
    sigma ~ dunif(0, 10)
  ), data=d
)

precis(m5.3)

plot(precis(m5.3))

# moving into some additional plotting methods now 

# building another model to plot residuals and other auxiliary data 

m5.4 <- map(
  alist(
    Marriage.s ~ dnorm(mu, sigma) , 
    mu <- a + b*MedianAgeMarriage.s , 
    a ~ dnorm(0, 10) , 
    b ~ dnorm(0, 1) , 
    sigma ~ dunif(0, 10)
  ), data=d
)

# To find residuals, use one predictor to model the other in a multivariate model 

# computing the expected value at MAP, for each State
mu <- coef(m5.4)['a'] + coef(m5.4)['b']*d$MedianAgeMarriage.s
# computing the residual for each state 
m.resid <- d$Marriage.s - mu

# and now plotting this 

plot(Marriage.s ~ MedianAgeMarriage.s, d, col=rangi2)
abline(m5.4)

# looping over the states 
for (i in 1:length(m.resid)) {
  x <- d$MedianAgeMarriage.s[i] # x location of line segment
  y <- d$Marriage.s[i] # observed endpoint of line segment
  # draw actual line segment
  lines(c(x,x), c(mu[i], y) , lwd=0.5, col=col.alpha("black", 0.7))
  
}

# preparing and plotting counterfactuals 

A.avg <- mean(d$MedianAgeMarriage.s)
R.seq <- seq(from=-3, to=3, length.out = 30)
pred.data <- data.frame(
  Marriage.s = R.seq, 
  MedianAgeMarriage.s=A.avg
)

# computing the counterfactual mean divorce 

mu <- link(m5.3, data=pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulating counterfactual divorce outcomes 

R.sim <- sim(m5.3, data=pred.data, n=1e4)

# and now plotting 

R.PI <- apply(R.sim, 2, PI)

# displaying predictions 
plot(Divorce ~ Marriage.s, data=d, type="n")
mtext("MedianAgeMarriage.s = 0")
lines(R.seq, mu.mean)
shade(mu.PI, R.seq)
shade(R.PI, R.seq)

# additional plotting
R.avg <- mean(d$Marriage.s)
A.seq <- seq(from=-3, to=3.5, length.out=30)
pred.data2 <- data.frame(
  Marriage.s = R.avg,
  MedianAgeMarriage.s=A.seq
)

mu <- link(m5.3, data=pred.data2)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

A.sim <- sim(m5.3, data=pred.data2, n=1e4)
A.PI <- apply(A.sim, 2, PI)

plot(Divorce ~ MedianAgeMarriage.s, data=d, type="n")
mtext("Marriage.s = 0")
lines(A.seq, mu.mean)
shade(mu.PI, A.seq)
shade(A.PI, A.seq)

?shade()

# posterior prediction plots, simulating predictions 

# call link without specifying new data, so it uses original data 
mu <- link(m5.3)

# summarizing samples across cases 
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulating observations , again using original data 
divorce.sim <- sim(m5.3, n=1e4)
divorce.PI <- apply(divorce.sim, 2, PI)

?PI

# now plotting things 
plot(mu.mean ~ d$Divorce, col=rangi2, ylim=range(mu.PI), xlab="Observed Divorce", ylab="Predicted Divorce")

abline(a=0, b=1, lty=2)
for(i in 1:nrow(d))
  lines(rep(d$Divorce[i], 2), c(mu.PI[1,i], mu.PI[2, i]) , col=rangi2)

identify(x=d$Divorce, y=mu.mean, labels=d$Loc, cex=0.8)

# code for some auxiliary plotting and analysis 

# compute residuals 

divorce.resid <- d$Divorce - mu.mean

# get ordering by divorce rate 

o <- order(divorce.resid)

# making the plot

dotchart(divorce.resid[o], labels=d$Loc[o], xlim=c(-6,5), cex=0.6)

# this is VERY nice and would be good for state-by-state analysis 

abline(v=0, col=col.alpha("black", 0.2))
for(i in 1:nrow(d)) {
  j <- o[i] # which State in order 
  lines(d$Divorce[j]-c(mu.PI[1,j], mu.PI[2,j]), rep(i,2))
  points(d$Divorce[j]-c(divorce.PI[1,j], divorce.PI[2,j]), rep(i,2), 
         pch=3, cex=0.6, col="gray")
}

# adding the third plot in the book's group 
N <- 100 # number of cases
x_real <- rnorm(N) # x_real as Gaussian with mean 0 and sd 1
x_spur <- rnorm(N, x_real) # x_spur as Gaussian with mean=x_real
y <- rnorm(N, x_real) # y as Gaussian with mean=x_real
d <- data.frame(y, x_real, x_spur) # binding all together in a data frame

head(d)

