#' @name splitLevel
#' @title create WOE using ctree to split level
#' @param x independent variable
#' @param y target variable
#' @param df data frame
#' @param minp minimum percentage in each level
#' @export
splitLevel <- function(formula, df, minp = 0.05)
{
  if (is.character(formula)) formula <- as.formula(formula)
  if (class(formula) != 'formula') stop("Formula is not correct.")
  if (any(!all.vars(formula) %in% names(df))) stop("Variables in formula is not in df.")

  elem <- as.list(formula)
  y <- eval(elem[[2]], df, parent.frame())
  x <- eval(elem[[3]], df, parent.frame())
  if (is.character(x)) formula <- as.formula(substitute(y ~ factor(x), list(y = elem[[2]], x = elem[[3]])))

  split <- ctree(formula, df, na.action = na.exclude, control = ctree_control(minbucket = minp * length(x)))
  bins <- width(split)
  #total one level with missing values
  if(bins == 1)
  {
    if(any(is.na(x)))
    {
    x_ <- is.na(x)
    band <- c('nonmissing', 'missing')
    }else{
      return(NULL)
    }
  }else{
    #deal with the ctree result and create `band` and split x into `x_`
    if(is.numeric(x))
    {
      breaks <- sort(unlist(sapply(1:length(split), function(i) split[[i]]$node$split$breaks)))
      x_ <- cut(x, breaks = unique(c(range(x, na.rm = T), breaks)), include.lowest = T, labels = F)
      band <- paste(c('', breaks), c(breaks, ''), sep = ' ~ ')
    }else if(is.ordered(x)){
      breaks <- sort(unlist(sapply(1:length(split), function(i) split[[i]]$node$split$breaks)))
      levels <- paste0("'", levels(x), "'")
      breaks <- cut(1:length(levels), breaks = unique(c(0, length(levels), breaks)), labels = F)
      band <- tapply(levels, breaks, paste0, collapse = ',')
      x_ <- breaks[match(as.character(x), levels(x))]
    }else{
      levels <- paste0("'", levels(x), "'")
      indexes <- unlist(sapply(1:length(split), function(i) split[[i]]$node$split$index))
      indexes[is.na(indexes)] <- ''
      indexes <- apply(matrix(indexes, ncol = length(levels), byrow = T), 2, paste0, collapse = '')
      x_ <- indexes[match(x, levels(x))]
      band <- tapply(levels, indexes, paste0, collapse = ',')
      #delete the levels in which the frequency equals zero
      #one reason is sampling which can reduce the frequency
      band <- band[names(band) != '']
    }
    if(any(is.na(x))) band <- c(band, 'missing')
  }

  #get the detail of split result
  splitResult <- table_matrix(x_, y, useNA = 'ifany')
  splitResult <- data.table(splitResult, keep.rownames = F)
  setnames(splitResult, c('CntGood', 'CntBad'))
  splitResult[, band := band]
  setcolorder(splitResult, c('band', 'CntGood', 'CntBad'))
  mode <- ifelse(is.numeric(x) || is.ordered(x), 'J', 'A')
  detail <- detail.woe(splitResult, mode)

  #linearity
  is.linear <- TRUE
  if(is.numeric(x) || is.ordered(x))
  {
    freqMatrix <- as.matrix(prop.table(table(x_, y, useNA = 'no'), 2))
    is.linear <- linearity(freqMatrix[,1], freqMatrix[,2]) < 1e-6
  }

  WoE_result <- list('summary' = data.frame('var'            = all.vars(elem[[3]]), #elem[[3]] maybe 'factor(x)'
                                            'class'          = class(x),
                                            'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                            'levels'         = nrow(splitResult) - any(is.na(x)),
                                            'IV'             = max(detail$IV),
                                            'IV_decrease'    = ifelse(class(x) %in% c("character", "factor"), NA, 0),
                                            'is.linear'      = is.linear,
                                            'is.suboptional' = FALSE,
                                            'method'         = 'ctree',
                                            'mode'           = mode,
                                            'detail'         = '',
                                            stringsAsFactors = F),
                     'detail' = detail,
                     'trace'  = NULL)
  class(WoE_result) <- 'woe.result'
  WoE_result
}
