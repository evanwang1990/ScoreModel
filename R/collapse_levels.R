freqMatrix <- function(x, y, levels)
{
  if(any(is.na(y))) stop("Target must not have NAs!\n")
  if(any(is.na(x)))
  {
    nas <- is.na(x)
    x <- x[!nas]
    y <- y[!nas]
  }

  if(is.numeric(x)) x <- cut(rank(x, ties.method = 'min'), breaks = levels, labels = F)

  freqMatrix <- matrix(table(x, y), ncol = length(unique(y)), dimnames = NULL)
  probMatrix <- apply(freqMatrix, 2, function(x) x/sum(x))
  return(list(freqMatrix = freqMatrix,
              probMatrix = probMatrix,
              is.ordered = is.numeric(x) || is.ordered(x),
              has.zero = any(freqMatrix == 0)))
}

getIV <- function(freqMatrix)
{
  if(freqMatrix$has.zero) stop("there are some zero cells!\n")

  probMatrix <- freqMatrix$probMatrix

  iv <- sum((probMatrix[,1] - probMatrix[,2]) * (log(probMatrix[,1]) - log(probMatrix[,2])))
  iv
}

getCstat <- function(freqMatrix)
{
  if(!freqMatrix$is.ordered) return(NULL)

  probMatrix <- freqMatrix$probMatrix
  c_stat <- cal_c_stat(probMatrix[,1], probMatrix[,2])
  c_stat
}

getXstat <- function(freqMatrix)
{
  if(!freqMatrix$is.ordered) return(NULL)

  probMatrix <- freqMatrix$probMatrix
  x_stat <- cal_x_stat(probMatrix[,1], probMatrix[,2])
  x_stat
}

collapseLevel <- function(x,
                          y,
                          levels,
                          method = c('iv', 'll'),
                          mode = c('A', 'J'),
                          do.trace = F)
{
  if(any(is.na(y))) stop("There are NAs in target variable!\n")
  if(any(is.na(x)))
  {
    nas <- is.na(x)
    x <- x[!nas]
    y <- y[!nas]
  }

  if(is.numeric(x)) x <- cut(rank(x, ties.method = 'min'), breaks = levels, labels = F)

  freqMatrix <- as.matrix(table(x, y))
  if(any(freqMatrix == 0)) stop("There are zero cells!\n")

  probMatrix <- apply(freqMatrix, 2, function(x) x/sum(x))



}
