#' Get_DEG
#'
#' Get differentail expression genes based on optional criteria
#'
#' @usage Get_DEG(x, criteria = c("p", "logfc"), threshold, is.sort = T)
#'
#' @param x a DE object containing results from DEA_test
#' @param criteria criteria to filter DEG, with options of p value or absolute value of log_FC
#' @param threshold threshold of the criteria to filter DEG
#' @param is.sort whether sorting genes based on criteria
#'
#' @return a DE object with DE_info containing interesting gene
#'
#'
#' @export
#'

Get_DEG = function(x, criteria, threshold, is.sort = T){
  res = x@DE_info
  if (criteria == "logfc") {
    index = abs(res$logfc) > threshold
    res = res[index, ]
    if (is.sort) res = res[order(abs(res$logfc), decreasing = T), ]
  } else if (criteria == "p"){
    index = res$p < threshold
    res = res[index, ]
    if (is.sort) res = res[order(res$p, decreasing = F), ]
  }
  x@DE_info = res
  return(x)
}
