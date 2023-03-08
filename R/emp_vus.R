#' @importFrom Rcpp evalCpp
#' @useDynLib emplikROCS, .registration = TRUE

#' @export
vus <- function(x, y, z, type = c("Ustat", "ties")){
  type <- match.arg(type)
  if (any(is.na(x)) | any(is.na(y)) | any(is.na(z))) return(NA)
  else {
    res <- switch(type,
                  Ustat = vusC_U(x, y, z),
                  ties = vusC_ties(x, y, z))
    return(res)
  }
}
