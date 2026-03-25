#### A fair coint toss

set.seed(25)
x <- sample(c("H", "T"), size = 20, replace = TRUE)

x

y <- table(x)
y

y_prop <- prop.table(table(x))
y_prop

barplot(y_prop)


#### Generating with different probabilities

y <- sample(0:2, size = 10000, replace = TRUE, prob=c(0.25, 0.60, 0.15))

head(y)

table(y)

prop.table(table(y))

barplot(table(y))

# * replace=TRUE to over ride the default sample without replacement
# * prob = to sample elements with different probabilities


#### Probability distributions in R

# The distributions of standard random variables are included in R.
# The main notation

# * r - random
# * q - quantile
# * d - density
# * p - distribution

# are combined with the distribution namesand specific
# parameters for these distributions

# * norm(mean, sd)
# * exp(rate)
# * pois(lambda)
# * t(df, ncp)
# * unif(min, max)
# * gamma(shape, scale)
# * binom(size, prob)
# * chisq(df, ncp)
# * logis(location, scale)

# to create functions' names in R.

# Note that normal distribution has two parameters - mean and variance -
# and is denoted by $X \sim N(\mu, \sigma^2)$, but in R one should write
# standard deviation instead of variance. Before working with distributions,
# check in the help parameters' values.
#
# See: Distributions in the stats package.
#
# Note: more distributions can be found in EnvStats, actuar and other packages.
#
# Thus for every distribution there are four commands. The commands for
# each distribution are prepended with a letter to indicate the
# functionality:
#
#   * d - returns values of the probability density function;

z <- seq(-3, 3, by=0.1)
w1 <- dnorm(x=z, 0, 1)
w1

plot(z, w1, type = 'l')
points(z, w1, col=2, cex=0.7)


# * p - returns the values of a distribution function
# F(x) = P(X \leq x)

w2 <- pnorm(q=z, 0, 1)
w2

plot(z, w2, type = 'l')
points(z, w2, col=2, cex=0.7)


# * q - returns the inverse cumulative density function (quantiles)
# p = F(x)
# x = F^{-1}(p)

u <- seq(0,1, by=0.05)
w3 <- qnorm(p=u, 0, 1)
w3

plot(u, w3, type = 'l')
points(u, w3, col=2, cex=0.7)

# * r - returns randomly generated numbers.

set.seed(222)
rn <- rnorm(n = 100, 0, 1)
rn

plot(1:100, rn)
abline(h=0, col="red")

hist(rn, freq = F)

library(rcompanion)
plotNormalHistogram(rn, prob = TRUE,
                     main = "Normal Distribution overlay on Histogram",
                     length = 100)

# Since in statistics we put a mathematical model for the real world situation, thus

# * population is defined by probability distribution
# * we need to estimate parameters to get all the information about the population
# * randomly generated numbers are simple random sample
# * in practice we have only one sample (the data set)
# * putting the mathematical model we intentionally do "mistakes",
# since approximate estimation
# is quicker to compute, than to derive a precise model for every situation.
# For example we model
# human's body mass or age using normal distribution, though body mass
# and height have only possitive
# values, while normal distribution is defined on real values.

#### Generate some processes

# (a) $T=10, \sigma=0.1$


T. <- 10
set.seed(20)
epsilon <- rnorm(T.,0,1)
i <- 1:T.
t <- i/T.
x <- (t)^2
plot(t,x, type="l")

sigma <- 0.1
y <- x+epsilon*sigma
plot(t, y, type="l")


# (b) $T=100, \sigma=0.1$

T. <- 100
set.seed(20)
epsilon <- rnorm(T.,0,1)
i <- 1:T.
t <- i/T.
x <- (t)^2
plot(t,x,type="l")

sigma <- 0.1
y <- x+epsilon*sigma
plot(t, y, type="l")


# (c) $T=10, \sigma=\sqrt{0.5}$

T. <- 10
set.seed(20)
epsilon <- rnorm(T.,0,1)
t <- sort(runif(T.,0,1))
sigma <- sqrt(0.5)
c. <- 5
x <- sin(c.*pi*t)
plot(x, type="l")

y <- x+epsilon*sigma
plot(t, y, type="l")

# (d) $T=1000, \sigma=\sqrt{0.5}$

T. <- 1000
set.seed(20)
epsilon <- rnorm(T.,0,1)
t <- sort(runif(T.,0,1))
sigma <- sqrt(0.5)
c <- 5
x <- sin(c*pi*t)
plot(x, type="l")

y <- x+epsilon*sigma
plot(t, y, type="l")

# (e) $T=1000, \sigma=0.2$

T. <- 1000
set.seed(20)
epsilon <- rnorm(T.,0,1)
t1 <- seq(0,1,len=T.)
t <- 2*t1
sigma <- 0.2
x <- sqrt(t)
plot(x, type="l")

y <- x+epsilon*sigma
plot(t, y, type="l")

# (f)

J <- 10
N <- 100
set.seed(50)
epsilon <- replicate(J, rnorm(N, 0, 1))
sigma <- 1
i <- 1:N
t <- i/100
#gammaj dimension 1XJ
gammaj <- replicate(J,rexp(1,1))
#xj dimension NxJ
xj <- replicate(J, t^2)
#in order to multiply transpose xj: 1xJ with JxN => 1xN,
#but transposed Nx1
x <- t(gammaj*t(xj))
y <- x+sigma*epsilon
matplot(t, y, type="l")

# (1)
t <- seq(0,2,length=101)

xt <- ifelse(t <= 1, (t-0.25)^2, 0.25*(t-1))

plot(t, xt, type = "b")

# (2)

set.seed(999)
t <- sort(runif(101,0,2))
xt <- ifelse(t <= 1, (t-0.25)^2, 0.25*(t-1))

plot(t, xt, type = "b")

#### Generate wiener process

n<-1000
T. <- 1
delta<-T./n
t <- seq(0,T.,delta)

set.seed(578)
epsilon <- rnorm(n,0,1)

W <- cumsum(c(0,epsilon*sqrt(delta)))

plot(t,W,type="l",
     xlab="t",ylab="W(t)",xlim=c(0,T.+delta))

B <- W - t * W[n+1]
plot(t,B,type="l",
     xlab="t",ylab="B(t)",xlim=c(0,T.+delta))


W_matr <- matrix(NA, ncol = 10, nrow = n+1)
B_matr <- matrix(NA, ncol = 10, nrow = n+1)
for(i in 1:10) {
  epsiloni <- rnorm(n,0,1)
  W_matr[,i] <- cumsum(c(0,epsiloni*sqrt(delta)))
  B_matr[,i] <- W_matr[,i] - t * W_matr[n+1, i]
}

matplot(W_matr, type="l")
matplot(B_matr, type="l")

library(pvar)
x <- rwiener(end = 1, frequency = 1000)
y <- rbridge(end = 1, frequency = 1000)

plot(x)
plot(y)

set.seed(555)
xt <- replicate(100, rwiener(end = 1, frequency = 1000))
matplot(xt, type = "l")



################################################################
library(refund)
data("DTI")
FA.cca <- DTI[complete.cases(DTI),]
FA.cca$ID <- factor(FA.cca$ID)
FA.cca

cca_df = refund.shiny::as_refundObj(FA.cca$cca)
library(ggplot2)
ggplot(cca_df, aes(x = index, y = value, group = id)) +
  geom_path(alpha = .5, color = "blue")



library(tidyr)
library(dplyr)
dti_subset <- DTI %>%
  filter(visit == 1) %>%
  drop_na(cca)

# Convert to long format
cca_long <- as.data.frame(dti_subset$cca) %>%
  mutate(id = dti_subset$ID, sex = dti_subset$sex, case = dti_subset$case,
         pasat = dti_subset$pasat) %>%
  pivot_longer(cols = starts_with("cca"),
               names_to = "position",
               values_to = "fa_value") %>%
  mutate(index = as.numeric(gsub("cca_", "", position))/93)

ggplot(cca_long, aes(x = index, y = fa_value, group = id, color = sex)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = sex), fun = mean, geom = "line", linewidth = 1.5) +
  labs(title = "CCA trackt profile by sex",
       x = "Index",
       y = "FA",
       color = "Sex") +
  theme_minimal()


ggplot(cca_long, aes(x = index, y = fa_value, group = id, color = case)) +
  geom_line() +
  labs(title = "CCA trackt profile by case",
       x = "Index",
       y = "FA",
       color = "Case") +
  theme_minimal()

library(plotly)
cca_na <- na.omit(cca_long)
cca_na$col <- floor(cca_na$pasat/10)
p <- plot_ly(cca_na,
        x = ~index,
        y = ~pasat,
        z = ~fa_value,
        split = ~id,
        #color = ~col,
        type = 'scatter3d',
        mode = 'lines',
        line = list(width = 2),
        opacity = 0.5) %>%
  layout(title = "FA profiles according to pasat",
         scene = list(xaxis = list(title = "distance along tract"),
                      yaxis = list(title = "PASAT"),
                      zaxis = list(title = "FA values")))

p

