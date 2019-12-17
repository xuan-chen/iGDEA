#' DE_test
#'
#' Perform differential expression analysis
#'
#' @usage DEA_test(x, group, option = c("Wilcox", "T", "NB"), verbose = TRUE)
#'
#' @param x a DE object from Normalize_data
#' @param group a colname from sample_info or a vetor of factor indexing two group for comparison
#' @param option statistical method to perform differential expression test
#' @param verbose Print the information of reference and process
#'
#' @return add a DE_info slot with gene information, log fold change and p value
#'
#' @export
#'

DEA_test = function(x, group, option = "Wilcox", verbose = T){
  dt = x@data$normalized

  if (length(group) == 1 & is.character(group)){
    index = as.factor(x@sample_info[,group])
    grp_name = group
  } else if (length(group) == nrow(x@sample_info)){
    index = as.factor(group)
    grp_name = "input_index"
  }

  ref = index[1]
  ref_index = which(index == ref)

  if (verbose) cat(paste0("----- Performing DEA with reference as ", grp_name, "=",  ref, " -----\n"))

  ##### Perform test #####
  dt1 = x@data$normalized[,ref_index]
  dt2 = x@data$normalized[,-ref_index]
  if (nrow(dt1) != nrow(dt2)){
    cat("Warning: two datasets should have same number of rows!\n")
    return(NULL)
  }

  x@DE_info = x@gene_info
  logfc = sapply(1:nrow(dt1), function(i) log(mean(dt1[i,]) / mean(dt2[i,])) )
  x@DE_info$logfc = logfc

  if (option == "Wilcox"){
    # print("wilcox")
    source("./R/wilcox_test.R")
    p = WilcoxTest(dt1, dt2)
  } else if (option == "T"){
    # print("t")
    source("./R/t_test.R")
    p = T_Test(dt1, dt2)
  } else if (option == "NB"){
    # setwd("./R")
    source("./R/main.R")
    countdata = x@data$normalized
    countdata = round(countdata)
    modelMatrix = cbind(rep(1, ncol(dt1)),
                        ifelse(index == ref, 0, 1))
    colnames(modelMatrix) = c("(Intercept)","Treatment")
    mobject = DEA_deseq(countdata,modelMatrix)
    res = data.frame(dispGeneEst = mobject$dispGeneEst, dispFit = mobject$dispFit,
                     dispMAP = mobject$dispersions, pval = mobject$Wald_Pvalue[,2])
    p = res$pval
  }

  x@DE_info$p = p

  if (verbose) cat("-----  Done! ----- \n")
  return(x)
}

