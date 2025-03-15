



# source("../MainFunctions/ebrahimi_estimator.r")


bootstrap_ebrahimi_estimator <- function(x, B){
  
  v.Bootstrap <- rep(0, B)
  
  for(b in 1:B) {
    same_values <- TRUE
    while (same_values) {
      bx <- sample(x, replace = TRUE)
      if (!all(bx == bx[1])) {
        same_values <- FALSE
        entropy_result <- ebrahimi_estimator(bx)
        if (is.finite(entropy_result)) {
          v.Bootstrap[b] <- entropy_result
          #cat("Entropy for Bootstrap Sample", b, ": ", v.Bootstrap[b], "\n")
        }
      }
    }
  }
  
  
  t <- ebrahimi_estimator(x)
  estimated_mean <- mean(v.Bootstrap[!is.na(v.Bootstrap)])
  
  return(2*t - estimated_mean)
  
}