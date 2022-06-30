#Code for Statistical Rethinking

library(rethinking)
#pg 52
p_grid <- seq(from=0, to=1, length.out=1000)
prob_p <- rep(1,1000)
prob_data <- dbinom(6, size=9, prob=p_grid)
posterior <- prob_data * prob_p
posterior <- posterior/sum(posterior)

samples <- sample(p_grid, prob=posterior, size = 1e4,replace = TRUE)
plot(samples)
dens(samples)

# add up posteriro probs where p<.5
sum(posterior[p_grid<.5])

quantile(samples,.8)

#pg 56
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1,1000)
likelihood <- dbinom(3, size=3, prob=p_grid)
posterior <- likelihood * prior
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, size=10000, replace = TRUE, prob=posterior)
PI(samples, prob=.5) #middle .5 of posterior
HPDI(samples, prob = .5) #densest part of posterior


#pg 58 (point estimation)
p_grid[which.max(posterior)]
chainmode(samples, adj = .01)
mean(samples)
median(samples)

#pg 60 (calculating loss function)
sum(posterior*abs(.5-p_grid)) #expected loss at .5
loss <- sapply(p_grid, function(d) sum(posterior * abs(d-p_grid))) #loss at all points
p_grid[which.min(loss)] #point of minimal loss


#pg 62 (likelihood of outcomes)
dbinom(0:2, size = 2, prob=.7)

dummy_w <- rbinom(100000, size=2, prob=.7) #simulate two tosses
table(dummy_w)/100000

dummy_w <- rbinom(100000, size=9, prob=.7)
simplehist(dummy_w, xlab="dummy water count")

#pg 66 (propagate uncertainty in parameters to predictions)
w <- rbinom(10000, size = 9, prob = samples)

#pg 72
pos <- replicate(1000, sum(runif(16,)))


#pg 74
growth <- replicate(10000, prod(1+ runif(12, 0, 0.1)))
dens(growth, norm.comp=TRUE)

big <- replicate(10000, prod(1+ runif(12,0,0.5)))
small <- replicate(10000, prod(1 + runif(12,0, 0.01)))

log.big <- replicate(10000, log(prod(1+runif(12, 0, 0.5))))




#pg. 78
w <- 6
n <- 9
p_grid <- seq(from=0, to=1, length.out = 100)
posterior <- dbinom(w,n,p_grid) * dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)


#pg. 79
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18,]

#pg. 82
curve(dnorm(x, 178, 20), from = 100, to = 250)
curve(dunif(x, 0,50), from = -10, to= 60)

sample_mu <- rnorm(10000, 178, 20)
sample_sigma <- runif(10000, 0, 50)
prior_h <- rnorm(10000, sample_mu, sample_sigma )
dens(prior_h)

  
#pg 85
mu.list <- seq( from = 150, 160, length.out = 100)
sigma.list <- seq(from=7,to=9, length.out = 100)
post <- expand.grid(mu=mu.list, sigma=sigma.list)
post$LL <- sapply(1:nrow(post), function(i)
            sum(dnorm(d2$height,post$mu[i], post$sigma[i], log = TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) + 
  dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))
contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)

sample.rows <- sample(1:nrow(post), size = 10000, replace = TRUE, prob = post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
plot(sample.mu, sample.sigma, cex=.5, pch=16, col = col.alpha(rangi2, 0.1))


#pg. 86
dens(sample.mu)
dens(sample.sigma)



#pg. 87

d3 <- sample(d2$height, size=20)








#pg87
mu.list <- seq(from=150, to=170, length.out=200) 
sigma.list <- seq(from=4, to=20, length.out=200)
post2 <- expand.grid(mu=mu.list, sigma=sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i)
  sum(dnorm(d3, mean=post2$mu[i], sd=post2$sigma[i], log=TRUE)))
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, TRUE) + dunif(post2$sigma, 0, 50, TRUE)
post2$prod <- exp(post2$prod - max(post2$prod))
sample2.rows <- sample(1:nrow(post2), size=1e4, replace=TRUE, prob=post2$prod)
sample2.mu <- post2$mu[sample2.rows]
sample2.sigma <- post2$sigma[sample2.rows]
plot(sample2.mu, sample2.sigma, cex=0.5, col=col.alpha(rangi2,0.1),
     xlab="mu", ylab="sigma", pch=16)

dens(sample2.sigma, norm.comp=TRUE)


dens(sample2.mu, norm.comp = TRUE)


#pg 88
flist <- alist(height ~ dnorm(mu,sigma),
               mu ~ dnorm(178,20),
               sigma ~ dunif(0,50))

m4.1 <- quap(flist, data = d2)
precis(m4.1)

m4.2 <- quap(
  alist(
    height ~dnorm(mu,sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0,50)
    
  ), data = d2
)
precis(m4.2)

#pg. 90
vcov(m4.1)
diag(vcov(m4.1))
cov2cor(vcov(m4.1))

post <- extract.samples(m4.1, n=100000)
head(post)


xbar <- mean(d2$weight)


m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm (0,1),
    sigma ~ dunif(0, 50))
  , data = d2)
  

mu <- link(m4.3)
str(mu)

weight.seq <- seq(from = 25, to = 70, by= 1)

mu <- link(m4.3, data=data.frame(weight=weight.seq))
str(mu)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

plot(height ~ weight, data = d2, xlim = c(24, 75) , col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)

#pg 108
sim.height <- sim (m4.3, data = list(weight=weight.seq))
str(sim.height)

height.PI <- apply(sim.height, 2, PI, prob=0.89)

mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.89)

plot(height ~ weight, d2, xlim = c(24, 75) ,col= col.alpha(rangi2, .5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)
shade (height.PI, weight.seq)


rm(list = ls())

data(Howell1)
d <- Howell1
plot(height ~ weight, d)

#pg 111
#quadratic
d$weight_s <- (d$weight - mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2
m4.5 <- quap(
        alist(
          height ~ dnorm(mu, sigma) ,
          mu <- a + b1*weight_s + b2*weight_s2 ,
          a ~ dnorm(178, 20),
          b1 ~ dlnorm(0,1) ,
          b2 ~ dnorm(0,1) , 
          sigma ~ dunif(0,50) )
        , data = d)

#pg 112
weight.seq <- seq(from = -2.2, to = 2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(m4.5, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=.89)
sim.height <- sim(m4.5, data = pred_dat)
height.PI <- apply(sim.height,2, PI, prob=.89)
        
plot(height ~ weight_s, d, col=col.alpha(rangi2, .5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

#pg113
#cubic
d$weight_s3 <- d$weight_s^3
m4.6 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,10),
    b3 ~ dnorm(0,10),
    sigma ~ dunif(0,50)), data=d)

weight.seq <- seq(from = -2.2, to = 2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2, 
                 weight_s3 = weight.seq^3)
mu <- link(m4.6, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=.89)
sim.height <- sim(m4.6, data = pred_dat)
height.PI <- apply(sim.height,2, PI, prob=.89)

plot(height ~ weight_s, d, col=col.alpha(rangi2, .5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

#pg 114
#natural scale for x-axis

plot(height ~ weight_s, d, col=col.alpha(rangi2, .5), xaxt="n")
at <- c(-2,-1,0,1,2)
labels <- at * sd(d$weight) + mean(d$weight)
axis(side = 1, at=at, labels = round(labels,1))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

rm(list = ls())
data(cherry_blossoms)
d <- cherry_blossoms
precis(d)


#pg 114 splines
plot(doy~year, data = d)

#pg 117
d2 <- d[complete.cases(d$doy),]
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0,1,length.out=num_knots))
library(splines)
B <- bs(d2$year,
        knots = knot_list[-c(1, num_knots)],
        degree = 3, intercept = TRUE)
plot(NULL, xlim=range(d2$year), ylim=c(0,1), xlab = "year", ylab="basis")
for(i in 1:ncol(B)) lines(d2$year, B[,i])

#pg. 119
m4.7 <- quap(
        alist(
          D ~ dnorm(mu, sigma),
          mu <- a + B %*% w , 
          a ~ dnorm(100,10), 
          w ~ dnorm(0,10),
          sigma ~ dexp(1)
        ), data = list(D=d2$doy, B=B) ,
        start = list(w = rep(0, ncol(B)))
)

post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)
plot(NULL, xlim = range(d2$year), ylim = c(-6,6),
     xlab = "year", ylab = "basis * weight")
for(i in 1:ncol(B))
  lines(d2$year, w[i] * B[,i])


mu <- link(m4.7)
mu_PI <- apply(mu, 2, PI, .97)
plot(d2$year, d2$doy, col=col.alpha(rangi2, .3), pch = 16)
shade(mu_PI, d2$year, col=col.alpha("black", .5))


#Chapter 5
#pg 125

rm(list = ls())
data("WaffleDivorce")
d <- WaffleDivorce
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

#Chapter 5
data(WaffleDivorce)
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

sd(d$MedianAgeMarriage)

m5.1 <- quap(
  alist(
    D~dnorm(mu,sigma),
    mu <- a +bA*A,
    a ~ dnorm(0,.2),
    bA~dnorm(0, .5),
    sigma~dexp(1)),
  data = d)
  
set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post=prior, data = list(A=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for(i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha("black", 0.4))

#pg 127
A_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.1, data = list(A=A_seq))
mu.mean <- apply( mu, 2, mean)
mu.PI <- apply(mu,2, PI)

plot(D ~ A, data = d, col = rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)

m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, .2),
    bM ~ dnorm(0, .5),
    sigma ~ dexp(1)),
  data = d)

M_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.2, data = list(M=M_seq))
mu.mean <- apply( mu, 2, mean)
mu.PI <- apply(mu,2, PI)

plot(D ~ M, data = d, col = rangi2, xlim = c(-2.5,2.5))
lines(M_seq, mu.mean, lwd=2)
shade(mu.PI, M_seq)

#pg 133
m5.3 <- quap(
  alist(
      D ~ dnorm(mu, sigma),
      mu <- a + bM * M + bA*A,
      a ~ dnorm(0, .2),
      bM ~dnorm(0, .5),
      bA ~dnorm(0, .5),
      sigma ~ dexp(1) ),
  data = d)
precis(m5.3)
plot(coeftab(m5.1, m5.2, m5.3))# , par=c("bA", "bM"))


#pg. 134
N <- 50
age <- rnorm(N)
mar <- rnorm(N, -age)
div <- rnorm(N, age)

m5.4 <- quap(
    alist(
      M ~ dnorm(mu, sigma), 
      mu <- a + bAM * A,
      a ~ dnorm(0, .2),
      bAM ~ dnorm(0, .5),
      sigma ~ dexp(1)),
    data = d)

#pg 136
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

#pg 138
mu <- link(m5.3)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

D_sim <- sim(m5.3, n=10000)
D_PI <- apply(D_sim, 2, PI)
plot(mu_mean ~ d$D, col=rangi2, ylim = range(mu_PI),
     xlab = "Observed divorce", ylab = "Predicted divorce")
abline(a=0,b=1,lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i], 2), mu_PI[,i], col=rangi2)

identify(x = d$D, y = mu_mean, labels = d$Loc)

#pg 141
rm(list = ls())

data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, .2),
    bM ~ dnorm(0, .5),
    bA ~ dnorm(0, .5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0, .2),
    bAM ~ dnorm(0, .5),
    sigma_M ~ dexp(1)),
      data = d)

A_seq <- seq(from = -2, to = 2, length.out = 30)

sim_dat <- data.frame(A=A_seq)
s <- sim(m5.3_A, data = sim_dat, vars=c("M", "D"))

#pg 142
plot(sim_dat$A, colMeans(s$D), ylim = c(-2,2), type = "l",
     xlab="manipulated A", ylab="counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")

sim2_dat <- data.frame(A = (c(20,30) - 26.1)/1.24)
s2 <- sim(m5.3_A, data = sim2_dat, vars=c("M", "D"))
mean(s2$D[,2] - s2$D[,1])

#pg 143
sim_dat <- data.frame(M=seq(from = -2, to = 2, length.out = 30), A = 0)
s <- sim(m5.3_A, data = sim_dat, vars = "D")
plot(sim_dat$M, colMeans(s), ylim=c(-2,2), type = "l",
     xlab = "manipulated M", ylab = "counterfactual D")
shade(apply(s,2,PI), sim_dat$M)
mtext("Total counterfactual effect of M on D")

#pg 144

library(rethinking)
data(milk)
d <- milk
str(d)

d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))


m5.5_draft <- quap(
      alist(
        K ~ dnorm(mu, sigma) , 
        mu <- a + bN*N,
        a ~ dnorm(0, 1),
        bN ~ dnorm(0,1),
        sigma ~ dexp(1)
      ), data = d
)

d$neocortex.perc

dcc <- d[complete.cases(d$K, d$N, d$M),]

m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma) , 
    mu <- a + bN*N,
    a ~ dnorm(0, 1),
    bN ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = dcc
)

prior <- extract.prior(m5.5_draft)
xseq <- c(-2,2)
mu <- link(m5.5_draft, post=prior, data = list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for(i in 1:50) lines(xseq, mu[i,], col=col.alpha("black", .3))

m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma) , 
    mu <- a + bN*N,
    a ~ dnorm(0, .2),
    bN ~ dnorm(0,.5),
    sigma ~ dexp(1)
  ), data = dcc
)

prior <- extract.prior(m5.5)
xseq <- c(-2,2)
mu <- link(m5.5, post=prior, data = list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for(i in 1:50) lines(xseq, mu[i,], col=col.alpha("black", .3))

#6.1, pg. 163
N <- 100
set.seed(909)
height <- rnorm(N,10,2)
leg_prop <- runif(N, 0.4, 0.5)
leg_left <- leg_prop*height + rnorm(N, 0, .02)
leg_right <- leg_prop*height + rnorm(N, 0, .02)
d <- data.frame(height, leg_left, leg_right)

m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl * leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dexp(1)
    ), data = d)
precis(m6.1)
plot(summary(m6.1))

post <- extract.samples(m6.1)
plot( bl ~ br, post, col=col.alpha(rangi2, .1), pch=16)

sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd = 2, xlab="sum of bl and br")

m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2,10),
    sigma ~ dexp(1))
    ,data = d
  )
precis(m6.2)



#pg 170
library(rethinking)
data(milk)
d <- milk
sim.coll <- function(r=.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat,
               sd=sqrt((1-r^2) * var(d$perc.fat)))
          m <- lm(kcal.per.g ~ perc.fat + x, data = d)
          sqrt(diag(vcov(m)))[2] #stdv of the parameter
}
rep.sim.coll <- function(r=.9, n=100){
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}
r.seq <- seq(from=0, to=.99, by = .01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n= 100))
plot(stddev~r.seq,type="l", col=rangi2,lwd=2, xlab="correlation")

#pg. 171
set.seed(71)
#number of plants
N <- 100

#simulate initial heights
h0 <- rnorm(N, 10, 2)

#assign treatments and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=.5 - treatment*.4)
h1 <- h0 + rnorm(N,5-3*fungus)
d <- data.frame(h0, h1, treatment, fungus)

#pg 172
sim_p <- rlnorm(10000, 0, .25)

m6.6 <- quap(
  
  alist(
    h1 ~dnorm(mu, sigma), 
    mu <- h0*p,
    p~dlnorm(0,.25),
    sigma ~dexp(1)
  ), data = d
)


m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 *p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, .2),
    bt ~ dnorm(0, .5),
    bf ~ dnorm(0, .5),
    sigma ~ dexp(1)
  ), data = d)

m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 *p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, .5),
    sigma ~ dexp(1)
  ), data = d)

##chapter 7
##pg. 194

sppnames <- c("afarensis", "africanus", "habilis", "boisei", "rudolfensis", 
              "ergaster", "sapiens")
brainvolcc <- c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5)
d <- data.frame(species=sppnames, brain=brainvolcc, mass = masskg)
d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

#pg. 196
m7.1 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b*mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0,1)
  ), data = d
)

set.seed(12)
s <- sim(m7.1)
r <- apply(s,2,mean) - d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1 - resid_var/outcome_var

R2_is_bad <- function(quap_fit){
  s <- sim(quap_fit, refesh = 0)
  r <- apply(s, 2, mean) - d$brain_std
  1 - var2(r)/var2(d$brain_std)
}

m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2,
    a~dnorm(0.5,1),
    b~dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ),
  data = d, start = list(b=rep(0,2)))

m7.3 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + 
      b[3]*mass_std^3,
    a~dnorm(0.5,1),
    b~dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ),
  data = d, start = list(b=rep(0,3)))


m7.4 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + 
      b[3]*mass_std^3 + b[4]*mass_std^4,
    a~dnorm(0.5,1),
    b~dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ),
  data = d, start = list(b=rep(0,4)))

m7.5 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + 
      b[3]*mass_std^3 + b[4]*mass_std^4 + 
      b[5]*mass_std^5,
    a~dnorm(0.5,1),
    b~dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ),
  data = d, start = list(b=rep(0,5)))

m7.6 <- quap(
  alist(
      brain_std ~ dnorm(mu, 0.001),
      mu <- a + b[1]*mass_std + b[2]*mass_std^2 + 
        b[3]*mass_std^3 + b[4]*mass_std^4  + 
        b[5]*mass_std^5 + b[6]*mass_std^6,
      a~dnorm(0.5,1),
      b~dnorm(0,10),
      log_sigma ~ dnorm(0,1)
    ),
    data = d, start = list(b=rep(0,6)))

#plots from pg. 200
{
post <- extract.samples(m7.1)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)
l <- link(m7.1, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)

post <- extract.samples(m7.2)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)
l <- link(m7.2, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)


post <- extract.samples(m7.3)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)
l <- link(m7.3, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)


post <- extract.samples(m7.4)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)
l <- link(m7.4, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)

post <- extract.samples(m7.5)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)
l <- link(m7.5, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)

post <- extract.samples(m7.6)
mass_seq <- seq(from=min(d$mass_std), to=max(d$mass_std),
                length.out=100)=
l <- link(m7.6, data=list(mass_std=mass_seq))
mu <- apply(l, 2, mean)
ci <- apply(l, 2, PI)
plot(brain_std ~ mass_std, data = d)
lines(mass_seq, mu)
shade(ci, mass_seq)

}

#pg 206
p <- c(.15, .7, .15)
-sum(p*log(p))

#pg. 226
set.seed(11)
WAIC(m6.7)

set.seed(77)
compare(m6.6, m6.7, m6.8, func=WAIC)
compare(m6.6, m6.7, m6.8, func=PSIS)

set.seed(91)
waic_m6.7 <- WAIC(m6.7, pointwise=TRUE)$WAIC
waic_m6.8 <- WAIC(m6.8, pointwise=TRUE)$WAIC
n <- length(waic_m6.7)
diff_m6.7_m6.8 <- waic_m6.7- waic_m6.8
sqrt(n*var(diff_m6.7_m6.8))
plot(compare(m6.6, m6.7, m6.8))

set.seed(92)
waic_m6.6 <- WAIC(m6.6, pointwise=TRUE)$WAIC
diff_m6.6_m6.8 <- waic_m6.6 - waic_m6.8
sqrt(n*var(diff_m6.6_m6.8))
set.seed(93)
compare(m6.6, m6.7, m6.8)@dSE


#pg 213 overthinking
#This takes A LONG TIME!
N <- 20
kseq <- 1:5
dev <- sapply(kseq, function(k){
  print(k);
  r <- replicate(10000, sim_train_test(N=N, k=k));
  c(mean(r[1,]), mean(r[2,]), sd(r[1,]), sd(r[2,]))
  
})



#pg 230

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)

m5.1 <- quap(
  alist(
    D~dnorm(mu,sigma),
    mu <- a +bA*A,
    a ~ dnorm(0,.2),
    bA~dnorm(0, .5),
    sigma~dexp(1)),
  data = d)


m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, .2),
    bM ~ dnorm(0, .5),
    sigma ~ dexp(1)),
  data = d)

m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M + bA*A,
    a ~ dnorm(0, .2),
    bM ~dnorm(0, .5),
    bA ~dnorm(0, .5),
    sigma ~ dexp(1) ),
  data = d)

set.seed(24071847)
compare(m5.1, m5.2, m5.3, func=PSIS)

set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3,pointwise = TRUE)
set.seed(24071847)
WAIC.m5.3 <- WAIC(m5.3, pointwise = TRUE)
plot(PSIS_m5.3$k, WAIC.m5.3$penalty, xlab = "PSIS Pareto k",
     ylab="WAIC penalty", col=rangi2, lwd=2)

m5.3t <- quap(
  alist(
    D ~ dstudent(2, mu, sigma),
    mu <- a + bM * M + bA*A,
    a ~ dnorm(0, .2),
    bM ~dnorm(0, .5),
    bA ~dnorm(0, .5),
    sigma ~ dexp(1) ),
  data = d)
  
##chapter 8
#pg. 242

library(rethinking)
data(rugged)
d <- rugged

#make log version of outcome
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)

#pg. 243
m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std -0.215),
    a ~ dnorm(1,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = dd)

set.seed(7)
prior <- extract.prior(m8.1)

#set up the plot dimensions

plot(NULL, xlim = c(0,1), ylim=c(0.5, 1.5),
     xlab="ruggedness", ylab="log GDP")
abline(h = min(dd$log_gdp_std), lty = 2)
abline(h = max(dd$log_gdp_std), lty = 2)

#draw 50 lines from the prior
rugged_seq <- seq(from=-.1, 1.1, length.out = 30)
mu <- link(m8.1, post=prior, data = data.frame(rugged_std=rugged_seq))
for (i in 1:50) lines(rugged_seq, mu[i,], col = col.alpha("black", .3))

sum(abs(prior$b) > .6)/length(prior$b)

#make more realistic priors
m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std -0.215),
    a ~ dnorm(1,.1),
    b ~ dnorm(0,.3),
    sigma ~ dexp(1)
  ), data = dd)
    
precis(m8.1)


#make a variable to index Africa (1) or not (2)
dd$cid <- ifelse(dd$cont_africa==1, 1, 2)

m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b*(rugged_std -0.215),
    a[cid] ~ dnorm(1,.1),
    b ~ dnorm(0,.3),
    sigma ~ dexp(1)
  ), data = dd)

compare(m8.1, m8.2)
precis(m8.2,depth=2)

post <- extract.samples(m8.2)
diff_a1_a2 <- post$a[,1] - post$a[,2]
PI(diff_a1_a2)

rugged.seq <- seq(from=-.1, to = 1.1, length.out = 30)
#computer mu over samples, fixing cid=2 and then cid=1
mu.NotAfrica <- link(m8.2,
                     data = data.frame(cid=2, rugged_std=rugged.seq))
mu.Africa <- link(m8.2, 
                  data= data.frame(cid=1, rugged_std=rugged.seq))
#summarize to means and intervals
mu.NotAfrica_mu <- apply(mu.NotAfrica,2, mean)
mu.NotAfrica_ci <- apply(mu.NotAfrica,2, PI, prob=.97)
mu.Africa_mu <- apply(mu.Africa,2, mean)
mu.Africa_ci <- apply(mu.Africa, 2, PI, prob=.97)

m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std -0.215),
    a[cid] ~ dnorm(1,.1),
    b[cid] ~ dnorm(0,.3),
    sigma ~ dexp(1)
  ), data = dd)

precis(m8.3, depth = 2)

compare(m8.1, m8.2, m8.3, func=PSIS)
plot(PSIS(m8.3, pointwise=TRUE)$k)

#plot Africa - cid=1
d.A1 <- dd[dd$cid==1,]
plot(d.A1$rugged_std, d.A1$log_gdp_std, pch=16, col=rangi2,
     xlab="ruggedness(standardized)", 
     ylab="log GDP (as proportion of mean)", 
     xlim = c(0,1))
mu <- link(m8.3, data = data.frame(cid=1, rugged_std=rugged_seq))
mu_mean <- apply(mu,2,mean)
mu_ci <- apply(mu,2, PI, prob=.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq,col=col.alpha(rangi2,.3))
mtext("African nations")


#plot Africa - cid=2
d.A0 <- dd[dd$cid==2,]
plot(d.A0$rugged_std, d.A0$log_gdp_std, pch=16, col=rangi2,
     xlab="ruggedness(standardized)", 
     ylab="log GDP (as proportion of mean)", 
     xlim = c(0,1))
mu <- link(m8.3, data = data.frame(cid=2, rugged_std=rugged_seq))
mu_mean <- apply(mu,2,mean)
mu_ci <- apply(mu,2, PI, prob=.97)
lines(rugged_seq, mu_mean, lwd=2)
shade(mu_ci, rugged_seq,col=col.alpha(rangi2,.3))
mtext("Non-African nations")

rugged_seq <- seq(from=-.2, to=1.2, length.out=30)
muA <- link(m8.3, data=data.frame(cid=1, rugged_std=rugged_seq))
muN <- link(m8.3, data=data.frame(cid=2, rugged_std=rugged_seq))
delta <- muA-muN



#fit a lm model to the data

lm.african_interaction <- lm(log_gdp ~  cont_africa + rugged_std + cont_africa*rugged_std, data = dd)
lm.african <- lm(log_gdp ~  cont_africa + rugged_std, data = dd)

summary(lm.african_interaction)
summary(lm.african)


#pg. 253

library(rethinking)
data(tulips)
d <- tulips
str(d)

d$blooms_std <- d$blooms/max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

a <- rnorm(10000, .5, 1)
sum(a < 0 | a > 1)/length(a)

a <- rnorm(10000, .5, .25)
sum(a < 0 | a > 1)/length(a)

m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a + bw*water_cent + bs*shade_cent,
    a ~ dnorm(.5, .25),
    bw ~ dnorm(0, .25),
    bs ~ dnorm(0, .25),
    sigma ~ dexp(1)
  ), data = d)

m8.5 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent,
    a ~ dnorm(.5, .25),
    bw ~ dnorm(0, .25),
    bs ~ dnorm(0, .25),
    bws ~ dnorm(0,.25), 
    sigma ~ dexp(1)
  ), data = d)

#pg. 258
par(mfrow=c(1,3)) 
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
       xlab="water", ylab = "blooms", pch=16, col=rangi2)
  mu <- link(m8.4, data= data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black", .3))
}

#same plot for m8.5
par(mfrow=c(1,3)) 
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
       xlab="water", ylab = "blooms", pch=16, col=rangi2)
  mu <- link(m8.5, data= data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black", .3))
}

#prior plots
set.seed(7)
prior <- extract.prior(m8.4)
par(mfrow=c(1,3)) 
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
       xlab="water", ylab = "blooms", pch=16, col=rangi2,
       main= paste("m8.4 prior: shade =", as.character(s)), 
       xaxt = "n", yaxt = "n")
  axis(1, at = c(-1,0,1))
  axis(2, at = c(0,.5,1))
  mu <- link(m8.4, data= data.frame(shade_cent=s, water_cent=-1:1), post=prior)
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black", .3))
}

#same plot for m8.5
set.seed(7)
prior <- extract.prior(m8.5)
par(mfrow=c(1,3)) 
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx], xlim=c(-1,1), ylim=c(0,1),
       xlab="water", ylab = "blooms", pch=16, col=rangi2)
  mu <- link(m8.5, data= data.frame(shade_cent=s, water_cent=-1:1), post=prior)
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black", .3))
}

#9.1
#

num_weeks <- 100000
positions <- rep(0,num_weeks)
current <- 10
for(i in 1:num_weeks){
  #record current position
  positions[i] <- current
  #flip coin to generate proposal
  proposal <- current + sample(c(-1,1), size = 1)
  #make sure he loops
  if(proposal < 1) proposal <- 10
  if(proposal > 10) proposal <- 1
  #move?
  prob_move <- proposal/current
  current <- ifelse(runif(1) < prob_move, proposal, current)
}

plot(1:100, positions[1:100])
plot(table(positions))

#pg 279

library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1,1,2)

#the same model as 8.3 but not using quap
#pg 280
data_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(data_slim)

m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1,.1),
    b[cid] ~ dnorm(0, .3),
    sigma ~ dexp(1)
  ), data = data_slim,chains=4, cores = 4
)

#chapter 10
#pg. 301

p <- list()
p$A <- c(0,0,10,0,0)
p$B <- c(0,1,8,1,0)
p$C <- c(0,2,6,2,0)
p$D <- c(1,2,4,2,1)
p$E <- c(2,2,2,2,2)

p_norm <- lapply(p, function(q) q/sum(q))
H <- sapply(p_norm,function(q) -sum(ifelse(q==0,0,q*log(q))))

ways <- c(1,90,1260,37800,113400)
logwayspp <- log(ways)/10

#pg 308

#candidate distributions
p <- list()
p[[1]] <- c(1/4, 1/4, 1/4, 1/4)
p[[2]] <- c(2/6, 1/6, 1/6, 2/6)
p[[3]] <- c(1/6, 2/6, 2/6, 1/6)
p[[4]] <- c(1/8, 4/8, 2/8, 1/8)

#expected value
sapply(p, function(p) sum(p*c(0,1,1,2)))

#entropy
sapply(p, function(p) -sum(p*log(p)))

#309
p <- 0.7
A <- c((1-p)^2 , p*(1-p), (1-p)*p, p^2)
-sum(A*log(A))

sim.p <- function(G=1.4){
    x123 <- runif(3)
    x4 <- ((G)*sum(x123)-x123[2]-x123[3])/(2-G)
    z <- sum(c(x123,x4))
    p <- c(x123, x4)/z
    list(H=-sum(p*log(p)), p=p)
}

H <- replicate(100000, sim.p(1.4))
dens(as.numeric(H[1,]), adj=.1)

entropies <- as.numeric(H[1,])
distributions <- H[2,]
max(entropies)
distributions[which.max(entropies)]

#chapter 11

library(rethinking)
data(chimpanzees)
d <- chimpanzees

d$treatment <- 1 + d$prosoc_left + 2*d$condition
xtabs(~ treatment + prosoc_left + condition, d)

m11.1 <- quap(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a,
    a ~ dnorm(0,10)),
  data=d)

set.seed(1999)
prior <- extract.prior(m11.1, n=10000)

p <- inv_logit(prior$a)
dens(p, adj = .1)

m11.2 <- quap(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a + b[treatment], 
    a ~dnorm(0, 1.5),
    b[treatment] ~ dnorm(0,10)
  ), data = d)
set.seed(1999)

prior <- extract.prior(m11.2, n=10000)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[,k]))

dens(abs(p[,1] - p[,2]), adj = .1)

m11.3 <- quap(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a + b[treatment], 
    a ~dnorm(0, 1.5),
    b[treatment] ~ dnorm(0,.5)
  ), data = d)
set.seed(1999)
prior <- extract.prior(m11.3, n=10000)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[,k]))
mean(abs(p[,1] - p[,2]))
dens(abs(p[,1] - p[,2]), adj = .1)

dat_list <- list(pulled_left = d$pulled_left,
                 actor = d$actor,
                 treatment = as.integer(d$treatment))
#pg 330
m11.4 <- ulam(
  alist(pulled_left ~ dbinom(1, p),
        logit(p) <- a[actor] + b[treatment],
        a[actor] ~ dnorm(0, 0.5),
        b[treatment] ~ dnorm(0, 0.5)),
  data = dat_list , chains= 4, cores = 4, log_lik = TRUE
)

precis(m11.4, depth=2)

post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)), xlim = c(0,1))

labs <- c("R/N", "L/N", "R/P", "L/P")
plot(precis(m11.4, depth=2, pars = "b"))


diffs <- list(db13 = post$b[,1] - post$b[,3],
              db24 = post$b[,2] - post$b[,4])
plot(precis(diffs))

pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[1,]

#pg 332
plot(NULL, xlim = c(1,28), ylim = c(0,1), xlab = "",
     ylab = "proportion left lever", xaxt = "n", yaxt = "n")
axis(2, at=c(0,0.5,1), labels = c(0,0.5,1))
abline(h=0.5, lty=2)
for (j in 1:7) abline(v=(j-1) * 4+4.5, lwd=.5)
for (j in 1:7) text((j-1)*4+2.5, 1.1, concat("actor", j), xpd=TRUE)
for (j in (1:7)[-2]){
  lines((j-1)*4+c(1,3),pl[j,c(1,3)], lwd = 2, col = rangi2)
  lines((j-1)*4+c(2,4), pl[j,c(2,4)], lwd=2, col=rangi2)
}
points(1:28, t(pl), pch=16, col= "white", cex=1.7)
points(1:28, t(pl), pch=c(1,1,16,16), col=rangi2, lwd=2)
yoff <- 0.01
text(1, pl[1,1]-yoff,"R/N", pos=1, cex=.8)
text(2, pl[1,2]+yoff,"L/N", pos=3, cex=.8)
text(3, pl[1,3]-yoff,"R/P", pos=1, cex=.8)
text(4, pl[1,4]+yoff,"L/P", pos=3, cex=.8)
mtext("observed proportions/n")

dat <- list(actor=rep(1:7, each = 4), treatment=rep(1:4,times=7))
p_post <- link(m11.4, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)

d$side <- d$prosoc_left + 1
d$cond <- d$condition + 1

dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  side = d$side,
  cond = d$cond
)

m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p) ,
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm(0,1.5),
    bs[side] ~ dnorm(0, .5),
    bc[cond] ~ dnorm(0,.5)
  ), data = dat_list2, chains = 4, cores = 4, log_lik = TRUE
)

compare(m11.5, m11.4, func=PSIS)

post <- extract.samples(m11.4, clean = FALSE)
str(post)

#Run STAN from direct STAN code
m11.4_stan_code <- stancode(m11.4)
m11.4_stan <- stan(model_code = m11.4_stan_code, data = dat_list, chains = 4, cores = 4)
compare(m11.4_stan, m11.4)

post <- extract.samples(m11.4)
mean(exp(post$b[,4]-post$b[,2]))

data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2* d$condition
d$side <- d$prosoc_left + 1
d$cond <- d$condition + 1

d_aggregated <- aggregate(
                  d$pulled_left,
                  list(treatment=d$treatment, actor = d$actor,
                       side = d$side, cond = d$cond),
                  FUN = sum)
          colnames(d_aggregated)[5] <- "left_pulls")

dat <- with(d_aggregated, list(
  left_pulls = left_pulls,
  treatment = treatment,
  actor = actor,
  side = side,
  cond = cond
))

m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom(18, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data = dat , chains = 4, cores = 4, log_lik = TRUE
)

compare(m11.6, m11.4, func = PSIS)

#pg 339
-2 * dbinom(6,9,0.2, log = TRUE)
-2*sum(dbern(c(1,1,1,1,1,1,0,0,0), 0.2, log = TRUE))


#340
library(rethinking)
data(UCBadmit)
d <- UCBadmit

dat_list <- list(
  admit = d$admit,
  applications = d$applications,
  gid = ifelse(d$applicant.gender == "male", 1, 2)
)

m11.7 <- ulam(
  alist(
    admit ~ dbinom(applications, p), 
    logit(p) <- a[gid],
    a[gid] ~ dnorm(0, 1.5)
  ), data = dat_list, chains = 4, cores = 4
)
precis(m11.7, depth = 2)

post <- extract.samples(m11.7)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a=diff_a, diff_p = diff_p))

#pg 342
postcheck(m11.7)
for(i in 1:6){
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines(c(x,x+1), c(y1,y2), col=rangi2, lwd=2)
  text(x+0.5, (y1+y2)/2 + 0.5, d$dept[x], cex = 0.8, col = rangi2)
  
}


dat_list$dept_id <- rep(1:6, each = 2)
m11.8 <- ulam(
  alist(
    admit ~ dbinom(applications, p), 
    logit(p) <- a[gid] + delta[dept_id],
    a[gid] ~ dnorm(0, 1.5),
    delta[dept_id] ~dnorm(0, 1.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 10000
)
precis(m11.8, depth = 2)

post <- extract.samples(m11.8)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a=diff_a, diff_p = diff_p))

pg <- with(dat_list, sapply(1:6, function(k) applications[dept_id==k]/sum(applications[dept_id==k])))
rownames(pg) <- c("male", "female")
colnames(pg) <- unique(d$dept)
round(pg, 2)

#pg 345


y <- rbinom(100000, 1000, 1/1000)
c(mean(y), var(y))

library(rethinking)
data(Kline)
d <- Kline
d

d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact == "high", 2, 1)

curve(dlnorm(x, 0, 10), from=0, to=100, n=200)
a <- rnorm(1000, 0, 10)
lambda <- exp(a)
mean(lambda)
curve(dlnorm(x, 3, .5), from =0, to=100, n=200)

N <- 100
a <- rnorm(N, 3, .5)
b <- rnorm(N, 0, 10)
plot(NULL, xlim=c(-2,2), ylim = c(0,100))
for(i in 1:N) curve(exp(a[i] + b[i] * x), add = TRUE, col=grau())

set.seed(10)
N <- 100
a <- rnorm(N, 3, .5)
b <- rnorm(N, 0, 0.2)
plot(NULL, xlim=c(-2,2), ylim = c(0,100))
for(i in 1:N) curve(exp(a[i] + b[i] * x), add = TRUE, col=grau())

x_seq <- seq(from = log(100), to = log(200000), length.out = 100)
lambda <- sapply(x_seq, function(x) exp(a+b*x))
plot(NULL, xlim = range(x_seq), ylim = c(0,500), xlab = "log population",
     ylab = "total tools")
for(i in 1:N) lines(x_seq, lambda[i,], col=grau(), lwd=1.5)

plot(NULL, xlim = range(exp(x_seq)), ylim = c(0,500), xlab = "population", 
     ylab = "total tools")
for (i in 1:N) lines(exp(x_seq), lambda[i,], col=grau(), lwd=1.5)


#pg 352

dat <- list(
  T = d$total_tools,
  P = d$P,
  cid = d$contact_id
)

#intercept only model
m11.9 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a,
    a~dnorm(3, 0.5)
  ), data = dat, chains = 4, cores = 4, log_lik=TRUE
)



#interaction model

m11.10 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a[cid] + b[cid]*P,
    a[cid]~dnorm(3, 0.5),
    b[cid]~dnorm(0, .2)
  ), data = dat, chains = 4, cores = 4, log_lik=TRUE
)

compare(m11.9, m11.10, func = PSIS)

k <- PSIS(m11.10, pointwise=TRUE)$k
plot(dat$P, dat$T, xlab="log population(std)", ylab="total tools",
     col=rangi2, pch=ifelse(dat$cid==1, 1, 16), lwd=2, ylim = c(0,75), cex=1+normalize(k))
ns <- 100
P_seq <- seq(from = -1.4, to = 3, length.out=ns)

#low contact (cid = 1)
lambda <- link(m11.10, data = data.frame(P=P_seq, cid=1))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

#high contact(cid = 2)
lambda <- link(m11.10, data = data.frame(P=P_seq, cid=2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty=1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

#Population on natural scale
plot(d$population, d$total_tools, xlab = "population",
     ylab = "total tools", col=rangi2, pch = ifelse(dat$cid==1, 1,16), lwd=2,
     ylim=c(0,75), cex=1+normalize(k))
ns <- 100
P_seq <- seq(from= -5, to = 3, length.out=ns)
pop_seq <- exp(P_seq*1.53 + 9)

lambda <- link(m11.10, data = data.frame(P=P_seq, cid=1))
lmu <- apply(lambda,2,mean)
lci <- apply(lambda, 2, PI)
lines(pop_seq, lmu, lty=2, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)
lambda <- link(m11.10, data = data.frame(P=P_seq, cid=2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(pop_seq, lmu, lty=1, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)

#overthinking box
dat2 <- list(T=d$total_tools, P=d$population, cid = d$contact_id)
m11.11 <- ulam(
  alist(
    T ~ dpois(lambda),
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1,1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data = dat2, chain = 4, cores = 4, log_lik=TRUE
)
#doesn't plot correctly
plot(d$population, d$total_tools, xlab = "population",
     ylab = "total tools", col=rangi2, pch = ifelse(dat$cid==1, 1,16), lwd=2,
     ylim=c(0,75), cex=1+normalize(k))
ns <- 100
P_seq <- seq(from= -5, to = 3, length.out=ns)
pop_seq <- exp(P_seq*1.53 + 9)

lambda <- link(m11.11, data = data.frame(P=P_seq, cid=1))
lmu <- apply(lambda,2,mean)
lci <- apply(lambda, 2, PI)
lines(pop_seq, lmu, lty=2, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)
lambda <- link(m11.11, data = data.frame(P=P_seq, cid=2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(pop_seq, lmu, lty=1, lwd=1.5)
shade(lci, pop_seq, xpd=TRUE)


#pg 357
num_days <- 30
y <- rpois(num_days,1.5)
num_weeks <- 4
y_new <- rpois(num_weeks,.5*7)

y_all <- c(y,y_new)
exposure <- c(rep(1,30), rep(7,4))
monastery <- c(rep(0,30), rep(1,4))
d <- data.frame(y = y_all, days = exposure, monastery = monastery)

d$log_days <- log(d$days)
m11.12 <- quap(
  alist(
    y~dpois(lambda),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1)
  ), data = d
)

post <- extract.samples(m11.12)
lambda_old <- exp(post$a)
lambda_new <- exp(post$a + post$b)
precis(data.frame(lambda_old, lambda_new))


#11.3.1 pg 360

N <- 500
income <- c(1,2,5)
score <- .5*income
p <- softmax(score[1], score[2], score[3])
career <- rep(NA,N) #this is an empty vector for choices
set.seed(34302)
for(i in 1:N) career[i] <- sample(1:3, size = 1, prob=p)

code_m11.13 <-"
data{
  int N; 
  int K; //number of possible careers
  int career[N]; //outcomes
  vector[K] career_income;
}
parameters{
vector[K-1] a; //intercepts
real<lower=0> b; //association of income with choice
}
model{
  vector[K] p;
  vector[K] s;
  a ~ normal(0, 1);
  b ~ normal(0, .5);
  s[1] = a[1] + b*career_income[1];
  s[2] = a[2] + b*career_income[2];
  s[3] = 0; //pivot
  p = softmax(s);
  career ~ categorical(p);

}
"

dat_list <- list(N=N, K=3, career=career, career_income = income)
m11.13 <- stan(model_code=code_m11.13, data=dat_list, chains=4, iter = 10000 )
precis(m11.13, 2)

post <- extract.samples(m11.13)

#logit scores
s1 <- with(post, a[,1] + b*income[1])
s2_orig <- with(post, a[,2] + b*income[2])
s2_new <- with(post, a[,2] + b*income[2]*2)

#compute probabilities for original and counterfactual 
p_orig <- sapply(1:length(post$b), function(i)
  softmax(c(s1[i], s2_orig[i], 0)))
p_new <- sapply(1:length(post$b), function(i)
  softmax(c(s1[i], s2_new[i], 0)))

p_diff <- p_new[2,] - p_orig[2,]
precis(p_diff)


N <- 500
family_income <- runif(N)
b <- c(-2, 0,2)
career <- rep(NA, N)
for (i in 1:N){
  score <- 0.5*(1:3) + b*family_income[i]
  p <- softmax(score[1], score[2], score[3])
  career[i] <- sample(1:3, size=1, prob=p) 
}

code_m11.14 <- "
  data{
    int N; // n observations
    int K; // n outcome values
    int career[N]; //outcome
    real family_income[N];
  }
  parameters{
    vector[K-1] a; //intercepts
    vector[K-1] b; //coefficients of family income
  }
  model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for(i in 1:N){
    for (j in 1:(K-1)) s[j] = a[j] + b[j]*family_income[i];
    s[K] = 0; //the pivot
    p = softmax(s);
    career[i] ~ categorical(p);
    }
  }
"

dat_list <- list(N=N, K=3, career=career, family_income = family_income)
m11.14 <- stan(model_code=code_m11.14, data = dat_list, chains = 4, iter= 10000)
precis(m11.14, 2)

#pg 364
library(rethinking)
data(UCBadmit)
d <- UCBadmit

m_binom <- quap(
  alist(
    admit ~ dbinom(applications,p),
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ), data = d
)

dat <- list(admit = d$admit, rej= d$reject)
m_pois <- ulam(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1,a2) ~ dnorm(0, 1.5)
  ), data = dat, chains = 3, cores = 3
)

inv_logit(coef(m_binom))

k <- coef(m_pois)
a1 <- k['a1']; a2 <- k['a2']
exp(a1)/(exp(a1) + exp(a2))