#' Partitionning the metacommunity CV2
#' 
#' Perform the partitioning of metacommunity CV2 presented in Segrestin & Leps (2021)
#' 
#' @usage cv2_decomp(data, nrand = NA)
#' 
#' @param data a data.frame. It must include one column named "Community" and 
#' one column named "Species". Values of annual biomass of individual 
#' populations must be organized in columns (one column per year). No extra 
#' column is allowed.
#' @param  nrand an integer (optional, default value is NA). Number of 
#' randomizations to estimate the range of pop.sync under the hypothesis of 
#' independent fluctuations between populations. When enabled, it increases the
#' computation time.

cv2_decomp <- function(data, nrand = NA){
  if(!"Species" %in% colnames(data)) stop("Column 'Species' not found in data")
  if(!"Community" %in% colnames(data)) stop("Column 'Community' not found in data")
  
  data$Species <- as.factor(data$Species)
  data$Community <- as.factor(data$Community)
  X <- data[, !colnames(data) %in% c("Species", "Community")]
  n <- sum(rowSums(X) == 0)
  if(n > 0) message(paste(n, "row(s) with no biomass were ommited"))
  data <- data[rowSums(X) != 0, ]
  X <- data[, !colnames(data) %in% c("Species", "Community")]
  
  n <- nrow(data)
  var_metacom <- var(colSums(X)) 
  mean2_metacom <- mean(colSums(X))^2
  CV2_metacom <- var_metacom / mean2_metacom
  
  mat_all <- cov(t(X)) 
  var_all <- diag(mat_all)
  pop_var <- sum(var_all) / mean2_metacom
  diag(mat_all) <- 0
  
  cov_all <- data.frame(com1 = rep(data$Community, n), sp1 = rep(data$Species, n),
                        com2 = rep(data$Community, each = n), sp2 = rep(data$Species, each = n),
                        cov = as.vector(mat_all) / mean2_metacom,
                        max = (rep(sqrt(var_all), n) * rep(sqrt(var_all), each = n)) / mean2_metacom)
  
  real_sp <- paste0(data$Species, data$Community)
  cov_all$in1 <- paste0(cov_all$sp1, cov_all$com2) %in% real_sp
  cov_all$in2 <- paste0(cov_all$sp2, cov_all$com1) %in% real_sp
  
  within <- cov_all[cov_all$sp1 == cov_all$sp2, ]
  within <- c(value = sum(within$cov),
              max = sum(within$max),
              Beta_MP = sum(within$max) - sum(within$cov),
              n = nrow(within))
  
  direct <- cov_all[cov_all$com1 == cov_all$com2 & cov_all$sp1 != cov_all$sp2, ]
  direct <- c(value = sum(direct$cov),
              max = sum(direct$max),
              delta = sum(direct$max) - sum(direct$cov),
              n = nrow(direct))
  
  indirect <- cov_all[cov_all$com1 != cov_all$com2 
                      & cov_all$sp1 != cov_all$sp2
                      & (cov_all$in1 | cov_all$in2), ]
  indirect <- c(value = sum(indirect$cov),
                max = sum(indirect$max),
                Beta_CCi = sum(indirect$max) - sum(indirect$cov),
                n = nrow(indirect))
  
  no <- cov_all[cov_all$com1 != cov_all$com2 
                & cov_all$sp1 != cov_all$sp2
                & !(cov_all$in1 | cov_all$in2),]
  no <- c(value = sum(no$cov),
          max = sum(no$max),
          Beta_CCno = sum(no$max) - sum(no$cov),
          n = nrow(no))
  
  res <- list(CV2 = CV2_metacom,
              pop_var = pop_var,
              Pop_sync_direct = direct,
              Pop_sync_intra = within,
              pop_sync_indirect = indirect,
              pop_sync_no = no)
  
  if (is.numeric(nrand)){
    rand <- matrix(NA, nrand, 4)
    for(i in 1:nrand){
      randomize <- apply(X, MARGIN = 1, sample, size = ncol(X))
      mat_all <- cov(randomize)
      diag(mat_all) <- 0
      
      cov_rand <- as.vector(mat_all) / mean2_metacom
      within_rand <- sum(cov_rand[cov_all$sp1 == cov_all$sp2])
      direct_rand <- sum(cov_rand[cov_all$com1 == cov_all$com2
                                  & cov_all$sp1 != cov_all$sp2])
      indirect_rand <- sum(cov_rand[cov_all$com1 != cov_all$com2 
                                    & cov_all$sp1 != cov_all$sp2
                                    & (cov_all$in1 | cov_all$in2)])
      no_rand <- sum(cov_rand[cov_all$com1 != cov_all$com2 
                              & cov_all$sp1 != cov_all$sp2
                              & !(cov_all$in1 | cov_all$in2)])
      rand[i, ] <- c(direct_rand, within_rand, indirect_rand, no_rand)
    }
    res <- c(res, rand = list(rand))
  }
  
  class(res) <- "cv.dec"
  return(res)                    
}

print.cv.dec <- function (x, ...) {
  cat("\nDecomposition of the metacommunity squared coefficient of variation")
  cat("\nSee Segrestin & Leps (2022)")
  cat("\n")
  pop_sync <- x$Pop_sync_intra[1] + x$Pop_sync_direct[1] + x$pop_sync_indirect[1] + x$pop_sync_no[1]
  if ("rand" %in% names(x)) {
    pop_sync_rand <- quantile(rowSums(x$rand), probs = c(0.025, 0.975))
    pop_sync_rand <- paste0("[", round(pop_sync_rand[1], 4), "; ", round(pop_sync_rand[2], 4), "]")
    cat(paste0("\nCV2 = ", round(x$CV2, 4), ", Pop.var = ", round(x$pop_var, 4), ", Pop.sync = ", round(pop_sync, 4), " ", pop_sync_rand))
  } else {
    cat(paste0("\nCV2 = ", round(x$CV2, 4), ", Pop.var = ", round(x$pop_var, 4), ", Pop.sync = ", round(pop_sync, 4)))
  }
  
  cat("\n")
  names_pop_sync <- c("Pop.sync[direct]", "Pop.sync[intra]", "pop.sync[indirect]", "pop.sync[no]")
  pop_sync_val <- round(unlist(lapply(x[3:6], "[", 1)), 4)
  names_hamm <- c("Delta", "Beta[MP]", "Beta[CCi]", "Beta[CCno]")
  pop_sync_hamm <- round(unlist(lapply(x[3:6], "[", 3)), 4)
  
  if ("rand" %in% names(x)){
    pop_sync_ind <- apply(x$rand, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    pop_sync_ind <- apply(pop_sync_ind, MARGIN = 2, 
                          function(x) paste0("[", round(x[1], 4), "; ", round(x[2], 4), "]"))
    df <- data.frame(paste0(names_pop_sync, " = ", pop_sync_val),
                     pop_sync_ind,
                     paste0(names_hamm, " = ", pop_sync_hamm))
    colnames(df) <- c("Segrestin & Leps (2022)", "Rand (95% CI)", "Hammond et al. (2020)")
  } else {
    df <- data.frame(paste0(names_pop_sync, " = ", pop_sync_val),
                     paste0(names_hamm, " = ", pop_sync_hamm))
    colnames(df) <- c("Segrestin & Leps (2022)", "Hammond et al. (2020)")
  }
  
  cat("\n")
  print(df, row.names = F, right = F)
  cat("\n")
}