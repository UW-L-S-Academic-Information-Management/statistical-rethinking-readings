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
