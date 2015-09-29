band.collapse <- function(x, x_, band_)
{
  if(!is.numeric(x)) return(band_)
  range_ <- sapply(band_, function(ranges) round(max(x[x_ %in% eval(parse(text = paste0(c('c(', ranges, ')'), collapse = '')))]), 2))
  res <- paste(c(' ', range_[-length(range_)]), c(range_[-length(range_)], ' '), sep = ' ~ ')
  res
}

band.split <- function()
{

}
