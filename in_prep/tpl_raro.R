tpl_change <- function(start, end, b){
  ts <- seq(from = start, to = end, by = sign(end-start))
  t_var <- c()
  for (t in 1:length(ts)) {
    t_var[t] <- rnorm(1, mean = ts[t], sd = sqrt(ts[t]^b))
  }
  plot(t_var)
  
  t_prev <- rnorm(length(ts), mean = start, sd = sqrt(start^b))
  t_ext <- rnorm(length(ts), mean = t_var[length(t_var)], sd = sqrt(t_var[length(t_var)]^b))
  tot <- c(t_prev, t_var, t_ext)
  
  plot(t_prev, col = 1,
       ylim = c(min(tot), max(tot)),
       xlim = c(0, length(tot)),
       main = paste("TPL b =", b),
       xlab = "time",
       ylab = "abundance")
  points(y = t_var, x = (length(t_prev)+1):((length(t_prev)+length(t_var))) , col = 2)
  points(y = t_ext, x = (((length(t_prev)+length(t_var)))+1):
           (((length(t_prev)+length(t_var)))+length(t_ext)),
         col = 3)
}
tpl_change(start = 100, end = 300, b = 1.2)
tpl_change(start = 300, end = 100, b = 1.2)


local_tpl <- function(x, ws) {
  slopes <- c()
  for (t in 1:(length(x)-ws) ) {
    var <- var(x[t:(t+ws-1)])
    mu <- mean(x[t:(t+ws-1)])
    slope <- log(var) / log(mu)
    res <- data.frame(var = var, mu = mu, slope = slope)
    slopes[[t]] <- res 
  }
  slopes <- do.call("rbind", slopes)
  return(slopes)
}



tpsim <- function(ny, pwr, mu, trend_mean, trend_sd){
  a <- rnorm(ny, mean = mu, sd = sqrt(mu^pwr))
  b <- c()
  trend <- rnorm(ny, mean = trend_mean, trend_sd)
  b[1] <- a[1]
  for (t in seq_along(a)) {
    b[t+1] <- b[t]*exp(trend[t])
  }
  
  slopes_a <- local_tpl(a, ws = 3)
  slopes_b <- local_tpl(b, ws = 3)
  
  par(mfrow = c(2,2))
  plot(slopes_a$slope, main = "tpl a")
  abline(h = pwr)
  plot(a, main = "a")
  plot(slopes_b$slope, main = "tpl b")
  abline(h = pwr)
  plot(b, main = "b")
  par(mfrow = c(1,1))
  c(mean(slopes_a$slope, na.rm = T), 
    mean(slopes_b$slope, na.rm = T))
}

tpsim(
  ny = 50,
  pwr = 1.5,
  mu = 100,
  trend_mean = 0,
  trend_sd = 0.1)

pp <- sim_mvcomm(n_sp = 30)
cc <- pp$sim_data
tpl(apply(cc, 2, var), apply(cc, 2, mean))
tt <- data.frame(ss = apply(cc, 2, var), 
                 mu = apply(cc, 2, mean))
bes <- apply(tt, 1, function(x){
  log(x[1]) / log(x[2])
})
mean(bes)
weighted.mean(bes, w = tt$mu/sum(tt$mu))


