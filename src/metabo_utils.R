#' Turn my long-format data.tables into a matrix fit for the KrumScripts
#'
#' @description The original KrumScript functions expect a matrix-like data.frame
#' where columns contain metabolites and rows contain the samples.
#' This function expects that there is a column named "sample" that uniquely
#' identifies a given metabolite per condition.
#' The actual metabolite name should be in the column specified via \code{which_ID}
#'
#' @param in.dt data.table in long format with a column for the value to be used
#' to fill the matrix \code{value} and columns for "sample","metabolite".
#' @param which_ID name for the column that contains the metabolite (or protein
#' or peptide) identifiers. Default: "metabolite"
#' @param which_value name of the column from which the values should be taken;
#' default: "value"
#'
#' @return data.frame with metabolites as columns and samples as rows
#'
#' @seealso \code{\link{imputeKNN}}. The \code{colnames} and \code{rownames} might
#' differ from the original ones as they will be run via \code{make.names}
make_KrumScript_matrix <- function(in.dt, which_ID = "metabolite", which_value = "value"){

  if(which_ID == which_value){stop("which_ID and which_value must specify different columns.")}

 # ABCutilities::check_columns(c("sample",which_ID, which_value), in.dt, "in.dt", "make_KrumScript_matrix")

  km <- data.table::dcast(in.dt[, c("sample",which_ID, which_value), with = FALSE],
                          as.formula(paste0("sample~", which_ID)),
                          value.var = which_value) %>% data.frame()
  rownames(km) <- make.names(km$sample)
  colnames(km) <- make.names(names(km))
  km$sample <- NULL
  return(km)
}



#' Adjust the meta info for metabolite data.table
#'
#' @description Very similar to \code{matrix_to_dt}, but more generally applicable
#' to metabolite data.tables.
#'
#' @param in.dt data.table with metabolite data that needs to be adjusted.
#' @param rm_Rname_cols Whether "Rnames" and "sample2", which are the results
#' of \code{make.names(metabolite)} and \code{make.names(sample)} should be
#' removed from the resulting \code{data.table}. Default: TRUE
#' \code{si_mets} object.
#'
#' @examples \dontrun{
#' vals.dt <- reshape2::melt(as.matrix(in.mat)) %>% data.table # get values
#' setnames(vals.dt, c("Var1","Var2", "value"), c("sample2", "Rnames", value_name))
#' }
#'
#' @seealso \code{\link{matrix_to_dt}}, \code{\link{imputeKNN_dt}}
#' @import data.table
#'
#' @export
#'
adjust_metainfo <- function(in.dt, rm_Rname_cols = TRUE){

  out.dt <- copy(in.dt)

  # make sure we have the column we need for the mergings
  if("sample2" %in% names(in.dt)){
    add_sample2 <- FALSE
  }else{
    if("sample" %in% names(in.dt)){
      add_sample2 <- TRUE
    }else{
      stop("You need to supply either a column named 'sample' or 'sample2' within in.dt.")
    }
  }

  if("Rnames" %in% names(in.dt)){
    add_Rnames <- FALSE
  }else{
    if("metabolite" %in% names(in.dt)){
      add_Rnames <- TRUE
    }else{
      stop("You need to supply either a column named 'metabolite' or 'Rnames' within in.dt.")
    }
  }

  if(add_sample2){
    out.dt$sample2 <- make.names(out.dt$sample)
  }
  if(add_Rnames){
    out.dt$Rnames <- make.names(out.dt$metabolite)
  }

  # get the original metabolite names back
  data("metabo_info", package = "BlenisMetData")
  vals.dt <- metabo_info[out.dt, on = "Rnames"]

  if(!add_Rnames){
    setnames(vals.dt, "ori.names", "metabolite")
  }else{
    setnames(vals.dt, "ori.names", "i.ori.names") # make sure it'll be removed
  }

  # get the values
  out.dt <- out.dt[vals.dt, on = c("sample2", "Rnames")]
  #out.dt <- out.dt[,-c("sample","Rnames","day","condition","replicate"), with = FALSE]

  # add the information about day, condition, replicates for metabolite-sample-
  # combinations that had originally been removed due to missing values
  data("si_mets", package = "BlenisMetData")
  out.dt <- si_mets[out.dt, on = "sample2"]

  # cosmetics to remove superfluous columns
  rm_cols <-  grep("^i\\.", names(out.dt), value = TRUE)
  if(rm_Rname_cols){
    rm_cols <- c("sample2", "Rnames", rm_cols)
  }
  out.dt <- out.dt[, -c(rm_cols), with = FALSE]

  return(out.dt)
}
