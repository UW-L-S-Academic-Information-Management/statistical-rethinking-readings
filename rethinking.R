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
  )
)