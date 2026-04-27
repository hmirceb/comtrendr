library(tidyverse)
library(broom)

sim_trend <- function(start, nyear, mu, sd) {
  n <- c(start)
  trends <- rnorm(nyear-1, mean = mu, sd = sd)
  for (t in seq_len(nyear)[-1]) {
    n[t] <- exp( log(n[t-1]) + trends[t-1] )
  }
  
  res <- list(sim_data = n,
              trend_mean = mean(trends),
              trend_sd = sd(trends)
              )
  
  return(res)
}


params <- expand.grid(n_sp = seq(4, 30, length.out = 7),
            trend_mean = seq(-0.10, 0.10, length.out = 7))

library(doFuture)
plan(multisession)
niter <- 1000
sims <- foreach(j = seq_len(nrow(params))) %dopar% {
  param_res <- foreach(i = seq_len(niter), .combine = "rbind") %do% {
    d <- sim_mvcomm(n_sp = params$n_sp[j], 
                    years = 1000, 
                    trend_mean = params$trend_mean[j],
                    trend_sd = 0,
                    corr = 0,
                    bimodal_trend = TRUE)
    res <- data.frame(v = sum(apply(d$sim_data, 2, var)),
                      v_sum = var(rowSums(d$sim_data)),
                      v2 = sum(apply(d$sim_data, 2, var_t2)),
                      v2_sum = var_t2(rowSums(d$sim_data)))
    res$iter <- i
    res$n_sp <- params$n_sp[j]
    res$trend_mean <- params$trend_mean[j]
    res
  }
}
plan(sequential)
sims <- do.call("rbind", sims)
sims$rmse_v <- sqrt((sims$v-sims$v_sum)^2)
sims$rmse_v2 <- sqrt((sims$v2-sims$v2_sum)^2)

save(sims, file = "in_prep/sims.RData")

# sims <- apply(params, 1, function(x){
#   param_res <- foreach(i = seq_len(niter), .combine = "rbind") %dopar% {
#     d <- sim_mvcomm(n_sp = x[[1]], 
#                     years = 1000, 
#                     trend_mean = x[[2]],
#                     trend_sd = 0,
#                     corr = 0,
#                     bimodal_trend = TRUE)
#     res <- data.frame(v = sum(apply(d$sim_data, 2, var)),
#                       v_sum = var(rowSums(d$sim_data)),
#                       v2 = sum(apply(d$sim_data, 2, var_t2)),
#                       v2_sum = var_t2(rowSums(d$sim_data)))
#     res$iter <- i
#     res$n_sp <- x[[1]]
#     res$trend_mean <- x[[2]]
#     res
#   }
#   }
#   )
# plan(sequential)
# sims <- do.call("rbind", sims)

sims %>% 
  pivot_longer(cols = v:v2_sum) %>%
  separate(name, sep = "_", into = c("method", "type")) %>%
  mutate(type = case_when(is.na(type) ~ "nosum",
                          .default = "sum")) %>%
  pivot_wider(values_from = value,
              names_from = type) %>% 
  ggplot()+
  geom_abline(slope = 1)+
  geom_point(aes(x = nosum, y = sum, color = method))+
  theme_classic()+
  facet_grid(n_sp~trend_mean, scales = "free")


library(patchwork)

(sims %>% 
  pivot_longer(cols = v:v2_sum) %>%
  separate(name, sep = "_", into = c("method", "type")) %>%
  mutate(type = case_when(is.na(type) ~ "nosum",
                          .default = "sum")) %>%
  pivot_wider(values_from = value,
              names_from = type) %>% 
  filter(method == "v") %>% 
  ggplot()+
  geom_abline(slope = 1)+
  geom_point(aes(x = nosum, y = sum), color = "blue")+
  theme_classic()+
  facet_grid(n_sp~trend_mean, scales = "free")+
    ggtitle("normal"))+
(sims %>% 
  pivot_longer(cols = v:v2_sum) %>%
  separate(name, sep = "_", into = c("method", "type")) %>%
  mutate(type = case_when(is.na(type) ~ "nosum",
                          .default = "sum")) %>%
  pivot_wider(values_from = value,
              names_from = type) %>% 
  filter(method == "v2") %>% 
  ggplot()+
  geom_abline(slope = 1)+
  geom_point(aes(x = nosum, y = sum), color = "red")+
  theme_classic()+
  facet_grid(n_sp~trend_mean, scales = "free")+
   ggtitle("ttqlv"))


sims %>% 
  ggplot()+
  geom_abline(slope = 1)+
  geom_point(aes(x = rmse_v, y = rmse_v2), alpha = 0.1)+
  theme_classic()+
  facet_grid(n_sp~trend_mean, scales = "free")

sims_anova <- sims %>% 
  pivot_longer(cols = rmse_v:rmse_v2)




errores <- sims %>% 
  mutate(diff = abs((v-v_sum)/v_sum),
         diff2 = abs((v2-v2_sum)/v2_sum),
         dev = (v-v_sum)^2 ,
         dev2 = (v2-v2_sum)^2 ) %>% 
  group_by(iter, n_sp, trend_mean) %>% 
  summarise(rmse = sqrt( sum(dev)/n()),
            rmse2 = sqrt( sum(dev2)/n()),
            mape = 100*sum(diff)/n(),
            mape2 = 100*sum(diff2)/n(),
            .groups = "drop") %>% 
  mutate(id = paste(sep = "_", n_sp, trend_mean)) 

rmse_anova <- errores %>% 
  pivot_longer(cols = c(rmse, rmse2))

errores %>% 
  ggplot(aes(x = rmse, y = rmse2))+
  geom_point()+geom_smooth(method = "lm")+
  facet_grid(n_sp ~ trend_mean, scales = "free")

errores %>% 
  ggplot(aes(x = mape, y =  mape2))+
  geom_point()+geom_smooth(method = "lm")+
  facet_grid(n_sp ~ trend_mean, scales = "free")

lm_rmse <- rmse_anova %>% 
  nest(data = c(value, name), .by = id) %>% 
  mutate(
    fit = map(data, ~lm(value ~ name, data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

by_id <- group_by(n_sp, trend_mean)
do(by_id,
   glance(
     lm(vale ~ name, data = rmse_anova)
   )
)


a <- sim_trend(start = 100, nyear =  10000, mu = 1/10000, sd = 0.01)
plot(a$sim_data)
cv_term(a$sim_data)
cv_term(a$sim_data, "two")
var(a$sim_data)
