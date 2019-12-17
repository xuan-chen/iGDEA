#' Create_DE_object
#'
#' Create a object in S4 class for differnetial expression analysis with
#'   expression matrix, sample information and gene annotation.
#'
#' @usage Create_DE_object(expr_mtx, sample_info = NULL, gene_info = NULL)
#'
#' @param expr_mtx expression matrix (n genes * p cells) for expression data
#' @param sample_info dataframe of p rows, storing sample information.
#'   If it is not specified, a data.frame is created with one
#'   column of colnames of expr_mtx as sample id.
#' @param gene_info dataframe of n rows, storing gene annotation information.
#'   If it is not specified, a data.frame is created with one
#'   column of rownames of expr_mtx as gene id.
#'
#' @return a S4 class containing all data and information
#'
#' @export
#'
#'
Create_DE_object = function(expr_mtx, sample_info = NULL, gene_info = NULL){
  e = as.matrix(expr_mtx)

  if (is.null(sample_info)) s = data.frame(sample_id = colnames(expr_mtx), stringsAsFactors = F) else
    s = data.frame(sample_id = colnames(expr_mtx), data.frame(sample_info), stringsAsFactors = F)
  if (is.null(gene_info)) g = data.frame(gene_id = rownames(expr_mtx), stringsAsFactors = F) else
    g = data.frame(gene_id = rownames(expr_mtx), data.frame(gene_info), stringsAsFactors = F)

  colnames(e) = rownames(s) = colnames(expr_mtx)
  rownames(e) = rownames(g) = rownames(expr_mtx)

  DE_object = setClass("DE Object", slots = list(data = "list", sample_info = "data.frame",
                                                 gene_info = "data.frame", DE_info = "data.frame"))
  new_object = DE_object(data = list(raw = e), sample_info = s, gene_info = g, DE_info = data.frame(NULL))

  return(new_object)
}
