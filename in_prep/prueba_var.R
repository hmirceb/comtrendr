n <- 20

m1 <- 20
d1 <- rnorm(n, mean = m1, sd = sqrt(m1))
d2 <- rnorm(n)

var(d1)
var(d2)
var(d1)+var(d2)
var(d1+d2)
(var(d1+d2)-var(d1)-var(d2))/2

Sigma <- matrix(c(1,0,0,1), 2, 2)
d <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)

var(d[,1])
var(d[,2])

cov(d[,1], d[,2])
(var(d[,1]+d[,2])-var(d[,1])-var(d[,2]))/2

var(d[,1]+d[,2])
detrending::var_t2(d[,1]+d[,2])
detrending::var_t2(d[,1])
detrending::var_t2(d[,2])
detrending::cov_term(d[,1], d[,2], term = "two")

plot(d[,1])
points(d[,2], col = "red")


d1 <- runif(n = n)
d2 <- runif(n = n)

var(d1)
var(d2)
cor(d1, d2)
