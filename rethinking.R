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


