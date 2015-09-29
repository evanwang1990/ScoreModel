#目前问题是：对于ordered factor应该使用breaks来处理，这个明天改进

splitLevel <- function(x,
                       y,
                       df,
                       minp = 0.05)
{
  if(is.character(x)) expr <- paste0(deparse(substitute(y)), ' ~ factor(', deparse(substitute(x)), ')')
  else expr <- paste(deparse(substitute(y)), '~', deparse(substitute(x)))

  split <- ctree(formula(expr), df, na.action = na.exclude, control = ctree_control(minbucket = minp * length(x)))
  bins <- width(split)
  #total one level including missing values
  if(bins == 1)
  {
    if(any(is.na(x)))
    {
    x_ <- is.na(x)
    band <- c('nonmissing', 'missing')
    }else{
      return(NULL)
    }
  }

  #deal with the ctree result and create `band` and split x into `x_`
  if(is.numeric(x))
  {
    breaks <- sort(unlist(sapply(1:length(split), function(i) split[[i]]$node$split$breaks)))
    x_ <- cut(x, breaks = c(range(x, na.rm = T), breaks), include.lowest = T, labels = F)
    band <- paste(c('', breaks), c(breaks, ''), sep = ' ~ ')
  }else{
    levels <- levels(x)
    indexes <- unlist(sapply(1:length(split), function(i) split[[i]]$node$split$index))
    indexes[is.na(indexes)] <- ''
    indexes <- apply(matrix(indexes, ncol = length(levels), byrow = T), 2, paste0, collapse = '')
    x_ <- as.character(x)
    for(indx in which(indexes != ''))
    {
      x_[x_ == levels[indx]] <- indexes[indx]
    }
    band <- sapply(sort(unique(indexes[indexes != ''])),
                   function(indx) paste0("'", paste(levels[which(indexes %in% indx)], collapse = "','"), "'"))
  }
  if(any(is.na(x))) band <- c(band, 'missing')

  #get the detail of split result
  splitResult <- dcast(as.data.frame(table(x_, y, useNA = 'ifany')), x_ ~ y, value.var = 'Freq', fill = 0)
  setDT(splitResult)
  setnames(splitResult, c('x_', 'CntGood', 'CntBad'))
  splitResult[, `:=`(x_     = NULL,
                     band   = band,
                     CntRec = CntGood + CntBad)
              ][, `:=`(PctRec   = paste0(round(CntRec / sum(CntRec) * 100, 2), '%'),
                       GoodRate = paste0(round(CntGood / CntRec * 100, 2), '%'),
                       BadRate  = paste0(round(CntBad / CntRec * 100, 2), '%'),
                       WoE      = ifelse(CntBad == 0, -4,
                                         ifelse(CntGood == 0, 4,
                                                round(log(CntBad / sum(CntBad)) - log(CntGood / sum(CntGood)), 4))))
                ][, `:=`(IV          = round(WoE * (CntBad / sum(CntBad) - CntGood / sum(CntGood)), 4),
                         WoE_barplot = barplot.woe(WoE))]
  setcolorder(splitResult, c('band', 'CntRec', 'PctRec', 'CntGood', 'CntBad', 'GoodRate', 'BadRate', 'WoE', 'IV', 'WoE_barplot'))
  if(!(is.numeric(x) || is.ordered(x))) setorder(splitResult, 'IV')
  splitResult <- rbind(splitResult,
                       data.frame(band        = 'Total',
                                  CntRec      = sum(splitResult$CntRec),
                                  PctRec      = '100%',
                                  CntGood     = sum(splitResult$CntGood),
                                  CntBad      = sum(splitResult$CntBad),
                                  GoodRate    = paste0(round(sum(splitResult$CntGood) / sum(splitResult$CntRec) * 100, 2), '%'),
                                  BadRate     = paste0(round(sum(splitResult$CntBad) / sum(splitResult$CntRec) * 100, 2), '%'),
                                  WoE         = 0,
                                  IV          = sum(splitResult$IV),
                                  WoE_barplot = ''))

  #linearity
  is.linear <- NA
  if(is.numeric(x) || is.ordered(x))
  {
    freqMatrix <- as.matrix(prop.table(table(x_, y, useNA = 'no'), 2))
    is.linear <- linearity(freqMatrix[,1], freqMatrix[,2]) < 1e-6
  }

  return(list('summary' = data.frame('var'             = deparse(substitute(x)),
                                     'class'          = class(x),
                                     'PctNA'          = round(sum(is.na(x)) / length(x), 3),
                                     'levels'         = nrow(splitResult) - 1 - any(is.na(x)),
                                     'IV'             = max(splitResult$IV),
                                     'IV_decrease'    = 0,
                                     'is.linear'      = is.linear,
                                     'is.suboptional' = F,
                                     'method'         = 'ct',
                                     'detail'         = ''),
              'detail'  = splitResult[1:(nrow(splitResult) - 1)]))
}
