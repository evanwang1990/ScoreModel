barplot.woe <- function(woe)
{
  woe <- scale(woe, center = T)
  woe <- round(woe / max(abs(woe)) * 5)
  res <- sapply(woe, function(x) ifelse(x < 0, paste0(c(rep(' ', 5 + x), rep('*', -x), '|', rep(' ', 5)), collapse = ''),
                                               paste0(c(rep(' ', 5), '|', rep('*', x), rep(' ', 5 - x)), collapse = '')))
  res
}
