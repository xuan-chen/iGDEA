#' Filter_Data
#'
#' Filter the genes with low expression and cells with low reads
#'
#' @usage Filter_Data(x, min_gene = 5, min_read = 5, verbose = T)
#'
#' @param x a DE object created from Create_simpleDE
#' @param min_gene filter out samples with expressed gene less than this
#' @param min_read filter out genes with total reads less than this
#' @param verbose Print the information of filtering
#'
#' @return return a filtered DE object with some genes and cells removed
#'
#'
#' @export
#'


Filter_Data = function(x, min_gene = 5, min_read = 5, verbose = T){

  mtx = x@data$raw

  s_index = colSums(mtx>0) > min_gene
  g_index = rowSums(mtx) > min_read

  s_filter = nrow(x@sample_info) - length(s_index)
  g_filter = nrow(x@gene_info) - length(g_index)

  x@data$raw = mtx[g_index, s_index]
  x@sample_info = x@sample_info[s_index, , drop = F]
  x@gene_info = x@gene_info[g_index, , drop = F]

  if (verbose){
    cat("\n")
    print(paste0(s_filter, " samples are filtered."))
    print(paste0(g_filter, " genes are filtered."))
  }

  return(x)
}
