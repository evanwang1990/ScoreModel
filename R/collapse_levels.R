#'
#'library(process)
#'df <- data.frame(x1 = sample(1:1000, 1e4, replace = T),
#'                 x2 = sample(c(letters[1:3], NA), 1e4, replace = T),
#'                 y  = c(sample(c(1, 1, 0), 5e3, replace = T), sample(c(1, 0, 0, 0), 5e3, replace = T)))
#'res <- bestCollapse(c('x1', 'x2'), y, df, 20, method = 'max_iv', mode = 'J', tracefile = 'trace.Rout', sqlfile = 'sql_code.sql')
#'collapseLevel(x = df$x2, y = df$y, levels = 20, method = 'max_iv', mode = 'A', minp = 0.05, sourcefile = 'test.R', sqlfile = 'test.sql')

#setcolorder
#对于其中的缺失值,求WOE
#WOE计算反了

band.collapse <- function(x, x_, band_)
{
  if(!is.numeric(x)) return(band_)
  range_ <- sapply(band_, function(ranges) round(max(x[x_ %in% eval(parse(text = paste0(c('c(', ranges, ')'), collapse = '')))]), 2))
  res <- paste(c(' ', range_[-length(range_)]), c(range_[-length(range_)], ' '), sep = ' ~ ')
  res
}

collapseLevel <- function(x, y, org.levels, method, mode, minp = 0.05, ...)
{
  if(is.character(x) || (is.factor(x) && !is.ordered(x)))
  {
    if(method == 'linear')
    {
      warning("Independent variable is character vector or unordered vector, 'linear' method is not appreciate\n The method is modified into 'max_iv' automatically.")
      method <- 'max_iv'
    }
    mode <- 'A'
  }else{
    mode <- 'J'
  }

  if(is.numeric(x))
  {
    x_ <- cut(rank(x, ties.method = 'min', na.last = "keep"), breaks = levels, labels = F)
  }else{
    x_ <- as.character(x)
  }

  freqMatrix <- as.matrix(table(x_, y))
  if(is.numeric(x))
  {
    labels <- rownames(freqMatrix)
  }else{
    labels <- paste0("'", rownames(freqMatrix), "'")
  }
  rownames(freqMatrix) <- labels

  # deal with the zero cells
  if(any(freqMatrix == 0))
  {
    freqMatrix_nonzero <- CollapseZeroCells(freqMatrix, matrix(NA, ncol = 2, nrow = nrow(freqMatrix)), mode = mode)
    # after processing, check if there is still zero cells
    if(any(freqMatrix_nonzero$freqMatrix == 0))
    {
      warning("There are zeros in some cells!\n")
      return(list("need_mo" = FALSE,
                  "main"    = data.frame(var         = deparse(substitute(x)),
                                         class       = class(x),
                                         NAs         = sum(is.na(x))/length(x),
                                         iv          = NA,
                                         levels      = NA,
                                         linear      = NA,
                                         suboptional = NA,
                                         method      = method,
                                         detail      = 'zero cells')))
    }

    #update labels of freqMatrix
    new_labels <- GetGroups(labels, freqMatrix_nonzero$label.trace[!is.na(freqMatrix_nonzero$label.trace[,1]),1], freqMatrix_nonzero$label.trace[!is.na(freqMatrix_nonzero$label.trace[,2]),2])
    freqMatrix <- freqMatrix_nonzero$freqMatrix
    rownames(freqMatrix) <- new_labels[new_labels != " "]
  }

  #check nrow
  nr <- nrow(freqMatrix)
  if(nr > 1)
  {
    trace <- matrix(nrow = nr,
                    ncol = 13,
                    dimnames = list(paste("Step", 0:(nr - 1)),
                                    c('Left', 'Right', 'minCount', 'max_iv', 'IV_decrease', 'X_stat', 'C_stat', 'Adjust_lift', 'Log_likehood', 'Prob(LR_Chi_Sq)', 'Z_score_of_log_odds_ratio', 'Prob(z_score)', 'Method')))
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 10] <- 1 - pchisq(-2 * (trace[-1, 9] - trace[1, 9]), 1:(nr - 1))
    trace[-1, 12] <- 1 - 2 * abs(pnorm(trace[-1, 11]) - 0.5)
    bin_iv <- trace[1, 12]
    trace[1, 12] <- 1
    best_indx <- max(min(which.max(trace[, 10] < 0.05) - 1, which.max(trace[, 12] < 0.05) - 1), which.max(trace[, 3] >= minp * length(x)))

    #get collapsed result----
    if(best_indx == 1){
      collapsed_result <- freqMatrix
    }else{
      collapsed_result <- combineResults(freqMatrix, trace[2:best_indx, 1], trace[2:best_indx, 2])
    }
    collapsed_result <- data.table(collapsed_result, keep.rownames = T)
    setnames(collapsed_result, c('band', 'CntGood', 'CntBad'))
    collapsed_result[, band := band.collapse(x, x_, band)]
    setorder(collapsed_result, c('band', 'CntGood', 'CntBad'))
    if (any(is.na(x)))
        collapsed_result <- cbind(collapsed_result,
                                  data.frame(band    = 'missing',
                                             CntGood = sum(is.na(x) & y == 0),
                                             CntBad  = sum(is.na(x) & y == 1)))
    collapsed_result <- detail.woe(collapsed_result, mode)

    args <- list(...)
    IV_ctree <- args[['IV_ctree']]
    WoE_result <- list('summary' = data.frame('var'             = deparse(substitute(x)),
                                              'class'          = class(x),
                                              'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                              'levels'         = nrow(collapsed_result) - 1 - any(is.na(x)),
                                              'IV'             = max(collapsed_result$IV),
                                              'IV_decrease'    = ifelse(is.null(IV_ctree), NA, round((max(collapsed_result$IV) - IV_ctree) / max(collapsed_result$IV), 3)),
                                              'is.linear'      = ifelse(mode == 'J', trace[best_indx, 8] < 1e-6, NA),
                                              'is.suboptional' = ifelse(mode == 'J', round((binary_IV - trace[nr - 1,4]) / bin_iv * 100, 1e-3), NA),
                                              'method'         = method,
                                              'mode'           = mode,
                                              'detail'         = ''),
                       'detail'  = collapsed_result[1:(nrow(collapsed_result) - 1)])
    #print trace----
    print.woe(trace = trace, best_index = best_indx, binary_IV = binary_IV, WoE_result = WoE_result)

  }else if(any(is.na(x))){  #the one level variables make sense only when they have missing values(just as two levels)
    stringcode <- stringCode(x, y, x_, deparse(substitute(x)), paste0("'",rownames(freqMatrix),"'"), !missing(sqlfile), method)
    writeLines(stringcode[1], sourcefile, sep = ',\n')
    if(!missing(sqlfile)) writeLines(stringcode[2], sqlfile)
    return(list("need_mo" = FALSE,
                "main"    = data.frame(var         = deparse(substitute(x)),
                                       class       = class(x),
                                       NAs         = sum(is.na(x))/length(x),
                                       iv          = 0,
                                       levels      = 1,
                                       linear      = ifelse(mode == 'J', trace[best_indx, 8] < 1e-6, NA),
                                       suboptional = 0,
                                       method      = method,
                                       detail      = 'one-level')))
  }
}


test <- function(...)
{
  list(...)
}
