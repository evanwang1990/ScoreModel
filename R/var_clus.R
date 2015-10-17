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

multipleR2 <- function(cluster, in_vars, all_vars, df)
{
  res <- matrix(NA, ncol = 4, nrow = length(in_vars), dimnames = list(in_vars, c('Cluster', 'SquaredR.own', 'SquaredR.next', '1-SquaredR Ratio')))
  res[, 1] <- cluster
  for (var in in_vars)
  {
    TSS <- sum(df[[var]]**2)

    if (length(in_vars) == 1)
      res[var, 2] <- 1
    else{
      formula1 <- as.formula(paste(var, '~', paste(setdiff(in_vars, var), sep = ' + ')))
      r1 <- residuals(lm(formula1, df))
      RSS1 <- sum(r1^2)
      res[var, 2] <- (TSS - RSS1) / TSS
    }

    if (length(in_vars) == length(var_names))
    {
      res[var, 3] <- 0
      res[var, 4] <- 0
    }else{
      formula2 <- as.formula(paste(var, '~', paste(setdiff(all_vars, in_vars), sep = ' + ')))
      r2 <- residuals(lm(formula2, df))
      RSS2 <- sum(r2^2)
      res[var, 3] <- (TSS - RSS2) / TSS
      res[var, 4] <- (1 - res[var, 2]) / (1 - res[var, 3])
    }
  }
  res
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
  clus[['node0']] <- list(vars     = 1:length(var_names),
                         need.deal = TRUE,
                         is.leaf   = TRUE)

  if (ncol(df) == 1) return(list(clusinfo = c(Cluster = 1, SquaredR.own = 1, SquaredR.next = 0, `1-SquaredR Ratio` = 0),
                                 var.exclude = list(var.char     = var.char,
                                                    var.nearzero = var.nearzero)))

  for(i in 1:maxiter)
  {
    deal.num <- 0
    for(j in 1:length(clus))
    {
      nodes <- clus[[j]]
      if (!nodes$need.deal) next

      twoclus <- twoCluster(df[,nodes$vars], maxeigen)

      if (is.null(twoclus) || all(twoclus))
      {
        clus[[j]][['need.deal']] <- FALSE
        next
      }

      clus[[paste0('node', i, '.1')]] <- list(vars      = nodes$vars[twoclus],
                                              need.deal = ifelse(sum(twoclus) == 1, FALSE, TRUE),
                                              is.leaf   = TRUE)
      clus[[paste0('node', i, '.2')]] <- list(vars      = nodes$vars[!twoclus],
                                              need.deal = ifelse(sum(!twoclus) == 1, FALSE, TRUE),
                                              is.leaf   = TRUE)
      clus[[j]][["need.deal"]] <- FALSE
      clus[[j]][['is.leaf']] <- FALSE
      deal.num <- 1
    }
    if (deal.num == 0) break
  }

  leaves <- which(unlist(lapply(clus, `[`, 3)))
  df <- as.data.frame(df)
  res <- lapply(1:length(leaves), function(i){node <- clus[[leaves[i]]]; multipleR2(i, var_names[node$vars], var_names, df)})
  res <- Reduce(rbind, res)
  return(list(clusinfo = res,
              var.exclude = list(var.char     = var.char,
                                 var.nearzero = var.nearzero)))
}
