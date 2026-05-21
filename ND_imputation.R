impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
  
  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")
  
  set.seed(seed)
  object <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })
  return(object)
}