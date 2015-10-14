band.collapse <- function(x, x_, band_)
{
  if(!is.numeric(x)) return(band_)
  range_ <- sapply(band_, function(ranges) round(max(x[x_ %in% eval(parse(text = paste0(c('c(', ranges, ')'), collapse = '')))]), 2))
  res <- paste(c(' ', range_[-length(range_)]), c(range_[-length(range_)], ' '), sep = ' ~ ')
  res
}

collapseLevel <- function(formula, df, org.levels = 20, method = c('max_iv', 'max_likehood', 'linear'), mode = 'J', minp = 0.05, ...)
{
  args <- list(...) #IV_ctree, skip.check
  IV_ctree <- args[['IV_ctree']]
  if (is.null(args[['skip.check']]))
  {
    if (is.character(formula)) formula <- as.formula(formula)
    if (class(formula) != 'formula') stop("Formula is not correct.")
    if (any(!all.vars(formula) %in% names(df))) stop("Variables in formula is not in df.")
  }

  elem <- as.list(formula)
  x <- eval(elem[[3]], df, parent.frame())
  y <- eval(elem[[2]], df, parent.frame())

  method <- match.arg(method)
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
    x_ <- cut(rank(x, ties.method = 'min', na.last = "keep"), breaks = org.levels, labels = F)
  }else{
    x_ <- as.character(x)
  }

  freqMatrix <- table_matrix(x_, y, useNA = 'no')
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
    freqMatrix <- CollapseZeroCells(freqMatrix, matrix(NA, ncol = 2, nrow = nrow(freqMatrix)), mode = mode)
    # after processing, check if there is still zero cells
    if(any(freqMatrix == 0))
    {
      freqMatrix <- data.table(freqMatrix, keep.rownames = T)
      setnames(freqMatrix, c('band', 'CntGood', 'CntBad'))
      freqMatrix[, band := band.collapse(x, x_, band)]
      setcolorder(freqMatrix, c('band', 'CntGood', 'CntBad'))
      if (any(is.na(x)))
        freqMatrix <- cbind(freqMatrix,
                                    data.frame(band    = 'missing',
                                               CntGood = sum(is.na(x) & y == 0),
                                               CntBad  = sum(is.na(x) & y == 1)))
      detail <- detail.woe(freqMatrix, mode)
      WoE_result <- list('summary' = data.frame('var'            = all.vars(elem[[3]]),
                                                'class'          = class(x),
                                                'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                                'levels'         = NA,
                                                'IV'             = NA,
                                                'IV_decrease'    = NA,
                                                'is.linear'      = NA,
                                                'is.suboptional' = NA,
                                                'method'         = method,
                                                'mode'           = mode,
                                                'detail'         = 'zero cells',
                                                stringsAsFactors = F),
                         'detail'  = detail,
                         'trace'   = NULL)
      class(WoE_result) <- 'woe.result'
      return(WoE_result)
    }
  }

  #check nrow
  nr <- nrow(freqMatrix)
  if(nr > 1)
  {
    row_names <- paste("Step", 0:(nr - 1))
    col_names <- c('Left', 'Right', 'minCount', 'IV', 'IV_decrease', 'X_stat', 'C_stat', 'Adjust_lift', 'Log_likehood', 'Prob(LR_Chi_Sq)', 'Z_score_of_log_odds_ratio', 'Prob(z_score)', 'Method')
    trace <- matrix(nrow = nr, ncol = 13, dimnames = list(row_names, col_names))
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 10] <- 1 - pchisq(-2 * (trace[-1, 9] - trace[1, 9]), 1:(nr - 1))
    trace[-1, 12] <- 1 - 2 * abs(pnorm(trace[-1, 11]) - 0.5)
    binary_IV <- trace[1, 12]
    trace[1, 12] <- 1
    best_indx <- max(min(which.max(trace[, 10] < 0.05) - 1, which.max(trace[, 12] < 0.05) - 1), which.max(trace[, 3] >= minp * length(x)))
    trace[1, 12] <- binary_IV #use trace[1, 12] to re-store binary_IV which will be useful in catLog function
    #star the best step in rownames
    row_names[best_indx] <- paste0(row_names[best_indx], '*')
    width <- max(nchar(row_names))
    row_names <- sapply(row_names, function(string) paste0(c(string, rep(' ', width - nchar(string))), collapse = ''))
    rownames(trace) <- row_names

    #get collapsed result----
    if(best_indx == 1){
      collapsed_result <- freqMatrix
    }else{
      collapsed_result <- combineResults(freqMatrix, trace[2:best_indx, 1], trace[2:best_indx, 2])
    }
    collapsed_result <- data.table(collapsed_result, keep.rownames = T)
    setnames(collapsed_result, c('band', 'CntGood', 'CntBad'))
    collapsed_result[, band := band.collapse(x, x_, band)]
    setcolorder(collapsed_result, c('band', 'CntGood', 'CntBad'))
    if (any(is.na(x)))
        collapsed_result <- rbind(collapsed_result,
                                  data.frame(band    = 'missing',
                                             CntGood = sum(is.na(x) & y == 0),
                                             CntBad  = sum(is.na(x) & y == 1)))
    detail <- detail.woe(collapsed_result, mode)

    WoE_result <- list('summary' = data.frame('var'            = all.vars(elem[[3]]),
                                              'class'          = class(x),
                                              'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                              'levels'         = nrow(detail) - 1 - any(is.na(x)),
                                              'IV'             = max(detail$IV),
                                              'IV_decrease'    = ifelse(is.null(IV_ctree), NA, round((max(detail$IV) - IV_ctree) / IV_ctree, 3)),
                                              'is.linear'      = ifelse(mode == 'J', trace[best_indx, 8] < 1e-6, TRUE),
                                              'is.suboptional' = ifelse(mode == 'J', round((binary_IV - trace[nr - 1,4]) / binary_IV * 100, 1e-3), NA),
                                              'method'         = method,
                                              'mode'           = mode,
                                              'detail'         = '',
                                              stringsAsFactors = F),
                       'detail'  = detail,
                       'trace'   = trace)
    class(WoE_result) <- 'woe.result'
    return(WoE_result)
  }else if(any(is.na(x))){
    #variable with two levels, missing and non-missing(after collapsing zeros)
    #if there is significant difference between missing and non-missing, output
    p.value <- chisq.test(table(is.na(x), y))$p.value
    if (p.value < 0.05)
    {
      x_ <- factor(is.na(x), levels = c(TRUE, FALSE), labels = c('missing', 'nomissing'))
      freqMatrix <- table_matrix(x_, y)
      freqMatrix <- data.table(freqMatrix, keep.rownames = T)
      setnames(freqMatrix, c('band', 'CntGood', 'CntBad'))
      detail <- detail.woe(freqMatrix, 'A')
      WoE_result <- list('summary' = data.frame('var'            = all.vars(elem[[3]]),
                                                'class'          = class(x),
                                                'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                                'levels'         = 1,
                                                'IV'             = max(detail$IV),
                                                'IV_decrease'    = ifelse(is.null(IV_ctree), NA, round((max(detail$IV) - IV_ctree) / IV_ctree, 3)),
                                                'is.linear'      = TRUE,
                                                'is.suboptional' = NA,
                                                'method'         = method,
                                                'mode'           = 'A',
                                                'detail'         = 'one level',
                                                stringsAsFactors = F),
                         'detail'  = detail,
                         'trace'   = NULL)
      class(WoE_result) <- 'woe.result'
      return(WoE_result)
    }
  }
}
