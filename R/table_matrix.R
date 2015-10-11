table_matrix <- function(x, y, useNA = c('ifany', 'no', 'always'))
{
  if (missing(useNA))
      useNA <- 'no'
  else
      useNA <- match.arg(useNA)
  table_ <- table(x, y, useNA = useNA)
  tableMatrix <- matrix(as.vector(table_), nrow = nrow(table_), dimnames = list(rownames(table_), colnames(table_)))
  tableMatrix
}
