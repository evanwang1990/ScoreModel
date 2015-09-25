barplot.woe <- function(woe)
{
  woe <- scale(woe, center = T, scale = F)
  woe[is.na(woe)] <- 0
  woe <- round(woe / max(abs(woe)) * 5)
  res <- sapply(woe, function(x) ifelse(x < 0, paste0(c(rep(' ', 5 + x), rep('*', -x), '|', rep(' ', 5)), collapse = ''),
                                               paste0(c(rep(' ', 5), '|', rep('*', x), rep(' ', 5 - x)), collapse = '')))
  res
}
