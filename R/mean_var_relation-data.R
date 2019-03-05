#' Relation between mean and sqrt(sd) of gene expression
#'
#'
#' @docType data
#'
#' @usage data(mean_var_relation)
#'
#' @format A list. The entries x and y are mean and sqrt(sd) of gene expression. 
#'
#' @keywords datasets
#'
#' @examples
#' data(mean_var_relation)
#' length(mean_var_relation)
#' names(mean_var_relation)
#' cbind(mean_var_relation$x, mean_var_relation$y)[1:2,]
"mean_var_relation"
