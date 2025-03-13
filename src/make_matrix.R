
#' Turn long data.table into matrix
#'
#' This is useful, e.g., for PCA and heatmaps.
#'
#' @param in.dt data.table in long format with a column for the value to be used
#' to fill the matrix (\code{which_value}) and columns for "sample" and
#' \code{which_ID}.
#' @param which_ID indicate the column with the variables for which the values
#' were observed, e.g. "metabolite", or "uniprot", or "tx_cluster"
#' @param log2_transform whether the data should be log2-transformed. Default:
#' FALSE. If set to TRUE, all "Inf" values will be replaced with NA.
#'
#' @details This returns a transposed version of \code{make_KrumScript_matrix}.
#'
#' @return data.frame with one column per sample and one row per instance of
#' \code{which_ID}, e.g., metabolites
make_matrix <- function(in.dt, which_ID = "metabolite", which_value = "value",
                        log2_transform = FALSE){

  if(which_ID == which_value){stop("which_ID and which_value must specify different columns.")}

  ABCutilities::check_columns(c("sample",which_ID, which_value), in.dt, "in.dt", "make_matrix")

  dt <- copy(in.dt)

  if(log2_transform){
    dt[, c(which_value) := log2(get(which_value))]
    dt[ get(which_value) == "-Inf", c(which_value) := NA]
      }

  km <- data.table::dcast(dt[, c("sample",which_ID, which_value), with = FALSE],
                          as.formula(paste0(which_ID, "~sample")),
                          value.var = which_value) %>% data.frame()
  colnames(km) <- make.names(names(km))
  rownames(km) <- make.names(km[, which_ID])
  km[which_ID] <- NULL

  return(km)
}
