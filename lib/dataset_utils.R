#dataset_utils

#' convert data frame to matrix
#' @param df data frame
#' @return a matrix
#' @description Converts a data frame to a matrix
#' @export
df.to.mat <- function(df) {
  row.col <- colnames(df)[1]
  df %>%
    as.data.frame() %>%
    select(1:3) %>%
    spread_(colnames(df)[2], colnames(df)[3]) %>%
    set_rownames(.[[row.col]]) %>%
    select_(paste("-", row.col)) %>%
    as.matrix()
}


#' convert matrix to data frame
#' @param mat matrix
#' @param row.name column label of resulting data frame from \code{mat} rows
#' @param col.name column label of resulting data frame from \code{mat} cols
#' @param dat.name column label of resulting data frame from \code{mat} data
#' @return a data frame
#' @description Converts a matrix to data frame
#' @export
mat.to.df <- function(mat, row.name="Row", col.name="Col", dat.name="Dat") {
  mat %>%
    as.data.frame %>%
    mutate_(.dots=c(~rownames(.)) %>% set_names(row.name)) %>%
    gather_(col.name, dat.name, colnames(mat)) %>%
    as.tbl
}

