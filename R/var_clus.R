twoCluster <- function(matrix, maxeigen)
{
  udv <- svd(matrix, nu = 0, nv = 2)
  pct <- udv$d / mean(udv$d)

  if (pct[2] < maxeigen) return(NULL)

  sd_y <- udv$d/sqrt(max(1, nrow(matrix) - 1))
  sd_x <- apply(matrix, 2, sd)
  cor_y1 <- sd_y[1] / sd_x * udv$v[,1]
  cor_y2 <- sd_y[2] / sd_x * udv$v[,2]

  return(abs(cor_y1) > abs(cor_y2))
}

varClus <- function(df, maxeigen, maxiter)
{
  var.char <- vector("character")
  var.nearzero <- vector("character")

  if (is.data.frame(df))
  {
    if(is.data.table(df)) setDF(df)
    var_class <- sapply(df, class)
    var.char <- names(var_class)[!var_class %in% c('numeric', 'integer')]
    var.int <- setdiff(names(var_class), var.char)
    if (length(var.int) == 0) stop("There are no Numeric variables in the data.")
    df <- as.matrix(df[, var.int])
  }

  if (is.character(df)) stop("The matrix should be numeric.")

  if (is.null(colnames(df))) stop("The matrix has no column names.")
  var.int <- colnames(df)

  df <- scale(df, scale = F)
  range_ <- apply(df, 2, function(x) max(x) - min(x))
  var.nearzero <- var.int[range_ < 10e-10]
  var_names <- setdiff(var.int, var.nearzero)
  if (length(var_names) == 0) stop("All the numeric variable are near zero")
  df <- df[,var_names]

  if (missing(maxiter)) maxiter <- length(df) - 1
  clus <- list()
  clus[['root']] <- list(vars = 1:length(var_names),
                         need.deal = 1)

  if (ncol(df) == 1) return(list(clus = clus,
                                 var.exclude = list(var.char     = var.char,
                                                    var.nearzero = var.nearzero)))

  for(i in 1:maxiter)
  {
    deal.num = 0
    for(j in length(clus))
    {
      nodes <- clus[[j]]
      if (nodes$need.deal == 0) next

      twoclus <- twoCluster(df[,nodes$vars], maxeigen)

      if (is.null(twoclus) || all(twoclus))
      {
        nodes$need.deal == 0
        next
      }

      clus[[paste0('node', i, '.1')]] <- list(vars      = nodes$vars[twoclus],
                                              need.deal = ifelse(sum(twoclus) == 1, 0, 1))
      clus[[paste0('node', i, '.2')]] <- list(vars      = nodes$vars[!twoclus],
                                              need.deal = ifelse(sum(!twoclus) == 1, 0, 1))
      clus[[j]][["need.deal"]] <- 0
      deal.num <- 1
    }
    if (deal.num == 0) break
  }
  return(list(clus = clus,
              var.exclude = list(var.char     = var.char,
                                 var.nearzero = var.nearzero)))
}
