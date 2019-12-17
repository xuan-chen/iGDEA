#' Differentail Expression Analysis using Student T-test
#'
#' Given two group of expression data, returns a vector of p values for each gene (row)
#'
#' @usage T_Test = function(data1, data2)
#'
#' @param data1 the reference matrix (m genes * n1 samples) for Wilcoxon Rank Sum Test
#' @param data2 the testing matrix (m genes * n2 samples) for Wilcoxon Rank SUm Test
#'
#' @details
#'
#' @export
#'
T_Test = function(data1, data2){
  if (nrow(data1) != nrow(data2)){
    cat("Warning: two datasets should have same number of rows!\n")
    return(NULL)
  }

  # data1 = ctr_dt; data2 = trt_dt;

  n1 = ncol(data1)
  n2 = ncol(data2)

  dt = cbind(data1, data2)

  res = apply(dt, 1, function(x){
    p = t.test(x[1:n1], x[(n1+1):(n1+n2)], "two.sided")$p.value
    return(p)
  })

  return(res)

}
