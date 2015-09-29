

注意：missing nonmissing的情况

stringCode <- function(x, y, x_, x_name, groups, SQLcode, method)
{
  #collapse groups
  groups <- paste0("c(", groups[groups != " "], ")")
  for(i in 1:length(groups))
  {
    x_[x_ %in% eval(parse(text = groups[i]))] <- i
  }

  #caculate range
  if(is.numeric(x))
  {
    range_ <- tapply(x, x_, max)
    if(SQLcode)
    {
      range_sql <- paste0(paste0(c('', range_[-length(range_)]), "<", x_name), "<=", c(range_[-length(range_)], ''))
      range_sql[1] <- substring(range_sql[1], 2)
      range_sql[length(range_sql)] <- substr(range_sql[length(range_sql)], 1, nchar(range_sql[length(range_sql)]) - 2)
    }
    range_ <- paste(c(-Inf, range_[-length(range_)]), c(range_[-length(range_)], Inf), sep = "<-")
  }else{
    range_ <- groups
    if(SQLcode) range_sql <- paste(x_name, 'in', substring(groups, 2))
  }

  #calculate WOE
  propMatrix <- as.matrix(prop.table(table(x_, y, useNA = 'ifany'), 2))
  woe <- round(log(propMatrix[,2] / propMatrix[,1]), unit = 1e-4)

  #string R code
  if(any(is.na(x))) range_ <- c(range_, "else")
  code <- paste0('"', paste(range_, woe, sep = "=", collapse = ";"), '"')
  code <- paste0(x_name, ifelse(method != 'mo', "_w = ", "_mw = "), paste0("recode(", paste(x_name, code, "as.factor.result = F", sep = ','), ")"))

  #string SQL code
  if(SQLcode)
  {
    sql_code <- paste('\twhen', range_sql, 'then', woe[!is.na(names(woe))], collapse = "\n")
    if(any(is.na(x))) sql_code <- paste(sql_code, '\n\telse', tail(woe,1))
    sql_code <- paste('case', sql_code, paste('as', paste0(x_name, ifelse(method == 'mo', '_mw', '_w')), ','), sep = '\n')
    return(c(code, sql_code))
  }
  return(code)
}
