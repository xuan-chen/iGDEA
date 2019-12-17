#' Normalize_Data
#'
#' Normalize raw data with optional normalization techniques
#'
#' @usage Normalize_Data(x, option = c("size_factor", "simplest", "median", "none"))
#'
#' @param x a DE object from Filter_Data
#' @param option option of normalization:
#'   "simplest": divided by mean reads of each cell
#'   "median": divided by median reads of each cell
#'   "size_factor": (default)divided by size factor, which is the median of the ratios of raw data and geometric mean
#'   "none": not performing normalization
#'
#' @return a DE object adding with a new element containing normalized data
#'
#'
#' @export
#'

Normalize_Data = function(x, option = "size_factor"){
  raw_data = x@data$raw

  if (option == "size_factor"){
    pseudo_ref = exp(rowMeans(log(raw_data + 1)))
    size_factor = apply(raw_data / pseudo_ref, 1, median)
    norm_data = raw_data / size_factor
  }else if (option == "simplest"){
    norm_data = raw_data / colMeans(raw_data)
  }else if (option == "median"){
    norm_data = raw_data / apply(raw_data, 2, median)
  }else if (option == "none"){
    norm_data = raw_data
  }

  x@data$normalized = norm_data

  return(x)
}
