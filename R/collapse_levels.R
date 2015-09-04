

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

  if(is.numeric(x)) x <- as.character(cut(rank(x, ties.method = 'min'), breaks = levels, labels = F))

  freqMatrix <- as.matrix(table(x, y))
  # check if there are zeros in cells
  if(any(freqMatrix == 0)) stop("There are zero cells!\n")


  #check nrow
  nr <- nrow(freqMatrix)
  if(nrow(freqMatrix) > 2)
  {
    trace <- matrix(nrow = nr - 1, ncol = 12)
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 9] <- 1 - pchisq(-2 * (trace[-1, 8] - trace[1, 9]), (nr - 2) : 1)
    trace[-1, 11] <- 1 - 2 * abs(pnorm(trace[-1, 10]) - 0.5)
    best_indx <- sum((trace[-1, 9] < 0.1) + (trace[-1, 11] < 0.1)  == 0) + 1

    #
    labels <- rownames(freqMatrix)
    groups <- GetGroups(labels, trace[2:best_indx, 1], trace[2:best_indx, 2])
  }







}
