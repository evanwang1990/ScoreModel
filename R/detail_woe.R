barplot.woe <- function(woe)
{
  woe <- scale(woe, center = T, scale = F)
  woe[is.na(woe)] <- 0
  woe <- round(woe / max(abs(woe)) * 5)
  res <- sapply(woe, function(x) ifelse(x < 0, paste0(c(rep(' ', 5 + x), rep('*', -x), '|', rep(' ', 5)), collapse = ''),
                                        paste0(c(rep(' ', 5), '|', rep('*', x), rep(' ', 5 - x)), collapse = '')))
  res
}

detail.woe <- function(woe_detail, mode)
{
  woe_detail[, `:=`(CntRec = CntGood + CntBad)
              ][, `:=`(PctRec   = paste0(round(CntRec / sum(CntRec) * 100, 2), '%'),
                       GoodRate = paste0(round(CntGood / CntRec * 100, 2), '%'),
                       BadRate  = paste0(round(CntBad / CntRec * 100, 2), '%'),
                       WoE      = ifelse(CntBad == 0, -4,
                                         ifelse(CntGood == 0, 4,
                                                round(log(CntBad / sum(CntBad)) - log(CntGood / sum(CntGood)), 4))))
                ][, `:=`(IV          = round(WoE * (CntBad / sum(CntBad) - CntGood / sum(CntGood)), 4),
                         WoE_barplot = barplot.woe(WoE))]
  setcolorder(woe_detail, c('band', 'CntRec', 'PctRec', 'CntGood', 'CntBad', 'GoodRate', 'BadRate', 'WoE', 'IV', 'WoE_barplot'))
  if(mode == 'A') setorder(woe_detail, 'IV')
  woe_detail <- rbind(woe_detail,
                      data.frame(band       = 'Total',
                                CntRec      = sum(woe_detail$CntRec),
                                PctRec      = '100%',
                                CntGood     = sum(woe_detail$CntGood),
                                CntBad      = sum(woe_detail$CntBad),
                                GoodRate    = paste0(round(sum(woe_detail$CntGood) / sum(woe_detail$CntRec) * 100, 2), '%'),
                                BadRate     = paste0(round(sum(woe_detail$CntBad) / sum(woe_detail$CntRec) * 100, 2), '%'),
                                WoE         = NA,
                                IV          = sum(woe_detail$IV),
                                WoE_barplot = ''))
  woe_detail
}
