tpl_change <- function(start, end, b){
  ts <- seq(from = start, to = end, by = sign(end-start))
  t_var <- c()
  for (t in 1:length(ts)) {
    t_var[t] <- rnorm(1, mean = ts[t], sd = sqrt(ts[t]^b))
  }
  plot(t_var)
  
  t_prev <- rnorm(151, mean = start, sd = sqrt(start^b))
  t_ext <- rnorm(151, mean = t_var[length(t_var)], sd = sqrt(t_var[length(t_var)]^b))
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

