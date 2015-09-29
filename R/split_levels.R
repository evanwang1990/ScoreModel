splitLevel <- function(x,
                       y,
                       df,
                       minp = 0.05)
{
  if(is.character(x)) expr <- paste0(deparse(substitute(y)), ' ~ factor(', deparse(substitute(x)), ')')
  else expr <- paste(deparse(substitute(y)), '~', deparse(substitute(x)))

  split <- ctree(formula, df, na.action = na.exclude, control = ctree_control(minbucket = minp * length(x)))
  bins <- width(split)
  #total one level including missing values
  if(bins == 1)
  {
    if(any(is.na(x)))
    {
      #----


    }else{
      return(NULL)
    }
  }

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
}
