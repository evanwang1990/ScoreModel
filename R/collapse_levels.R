collapseNames <- function(names, i, j)
{
  names[i] <- paste0(names[c(i, j)], collapse = '+')
  names[-j]
}

collapseLevel <- function(x,                                # independent variable
                          y,                                # target
                          levels,                           # initial levels
                          method = c('iv', 'll', 'mo'),     # collapse methods
                          mode = 'J',                       # collapse mode
                          do.trace = T)                     # need output collapsing step?
{
  # deal with NAs
  if(any(is.na(y))) stop("There are NAs in target variable!\n")
  if(any(is.na(x)))
  {
    nas <- is.na(x)
    x <- x[!nas]
    y <- y[!nas]
  }

  # deal with `method` and `mode`
  method <- match.arg(method)
  if(is.character(x) || (is.factor(x) && !is.ordered(x)))
  {
    if(method == 'mo')
    {
      warning("Independent variable is character vector or unordered vector, 'mo' method is not appreciate\n The method is modified into 'iv' automatically.")
      method <- 'iv'
    }
    mode <- 'A'
  }

  if(is.numeric(x)) x <- cut(rank(x, ties.method = 'min'), breaks = levels, labels = F)

  freqMatrix <- as.matrix(table(x, y))
  # check if there are zeros in cells
  if(any(freqMatrix == 0)) stop("There are zero cells!\n")


  #check nrow
  if(nrow(freqMatrix) > 2)
  {
    res <- collapse(freqMatrix, method, mode)
  }







}
