n <- 20
pwr <- 1.8
d <- seq_len(n)+n
for (t in seq_along(d)) {
  d[t] <- rnorm(1, mean = d[t], sd = sqrt(d[t]^pwr) )
}
plot(d)
abline(h = mean(d))

d2 <- seq_len(n)+n
for (t in seq_along(d2)) {
  d2[t] <- rnorm(1, mean = d2[t], sd = sqrt(d2[t]^pwr) )
}

d3 <- c(d, d2)
plot(d3)
abline(h = mean(d3))
abline(h = mean(d2), col = "red")
abline(h = mean(d), col = "blue")
