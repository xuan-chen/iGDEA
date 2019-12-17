#' Differnetial Expression Analysis using Wilcoxon Rank Sum Test
#'
#' Given two group of expression data, returns a vector of p values for each gene (row)
#'
#' @usage WilcoxTest = function(data1, data2)
#'
#' @param data1 the reference matrix (m genes * n1 samples) for Wilcoxon Rank Sum Test
#' @param data2 the testing matrix (m genes * n2 samples) for Wilcoxon Rank SUm Test
#'
#' @details
#'
#' @export

WilcoxTest = function(data1, data2){
  if (nrow(data1) != nrow(data2)){
    cat("Warning: two datasets should have same number of rows!\n")
    return(NULL)
  }

  # data1 = ctr_dt; data2 = trt_dt;

  n1 = ncol(data1)
  n2 = ncol(data2)

  dt = cbind(data1, data2)

  res = apply(dt, 1, function(x){
    p = WilcoxTest_Vec(x[1:n1], x[(n1+1):(n1+n2)])
  })

  return(res)

}

#' Wilcox rank sum test on two vector
#'
#' Given two vector, returns p value for Wilcox rank sum test
#'
#' @usage WilcoxTest_Vec(x, y)
#'
#' @param x reference vector for Wilcox rank sum test
#' @param y testing vector for Wilcox rank sum test
#'
#' @export
#'
WilcoxTest_Vec = function(x, y){

  n1 = length(x)
  n2 = length(y)
  n = n1 + n2
  r = rank(c(x,y))
  U = sum(r[1:n1]) - n1 * (n1 + 1) / 2
  # print(U)

  APX_threshold = 0
  APX = (n1 >= APX_threshold) && (n2 >= APX_threshold) # if approximated by normal distribution

  if (APX){
    # Approximated by Normal Distribution
    TIE = (length(unique(r)) != length(r)) # if any ties rank

    if (!TIE){
      # no tie
      sigma = sqrt(n1 * n2 * (n + 1) / 12)
    } else{
      # with tie
      TIE_table = table(r)
      sigma = sqrt(n1 * n2 / 12 * (n+1 - sum(TIE_table ^ 3 - TIE_table) / n / (n-1)) )
    }

    Z = (U - n1 * n2 / 2) / sigma
    p = 2 * min(pnorm(Z), pnorm(-Z))

  } else{
    # Testing using Wilcox Distribution
    p = pwilcox(U, n1, n2)
    p = 2 * min(p, 1-p)
  }

  return(p)
}
