#注意：missing nonmissing的情况
#注意：变量长度太长时，需要截断加密，制作对应表map
str_format <- function(str, ...)
{
  args <- as.list(unlist(list(...), recursive = T))
  if (length(args) == 0) return(str)
  if (is.null(names(args))) names(args) <- 1:length(args)
  for (i in 1:length(args))
  {
    str <- gsub(paste0('#', names(args[i])), args[[i]], str)
  }
  str
}

`%+%` <- function(x, y)
{
  paste0(x, y)
}


createConditions <- function(varname, class, band, type)
{
  if ('nonmissing' %in% band)
  {
    if (type == 'sql') condtions <- paste(varname, "is", ifelse(band == 'missing', '', 'not'), 'NULL')
    else conditions <- paste(ifelse(band == 'missing', '', '~'), "missing(", varname, ")")
  }else if (class %in% c("numeric", "integer")){
    breaks <- sapply(band[band != 'missing'], function(str) strsplit(str, ' ~ ', fixed = T)[[1]][1])
    breaks <- breaks[!is.na(breaks)]
    conditions <- paste(breaks, '<', varname, '<=', c(breaks[-1], ' '))
    conditions[1] <- substring(conditions[1],4)
    conditions_last <- conditions[length(condtions)]
    conditions[length(conditions)] <- substr(conditions_last, 1, nchar(conditions_last) - 5)
    if ('missing' %in% band) conditions <- c(conditions, ifelse(type == 'sql', paste(varname, "is NULL"), paste('missing(', varname, ')')))
  }else{
    conditions <- paste(varname, "in (", band[band != "missing"], ")")
    if ("missing" %in% band) conditions <- c(conditions, ifelse(type == 'sql', paste(varname, "is NULL"), paste('missing(', varname, ')')))
  }
  conditions
}

toSql <- function(varname, class, band, WoE)
{
  conditions <- createConditions(varname, class, band, 'sql')
  sql <- paste0(paste(c("case when", rep("     when", length(conditions) - 1)), conditions, "then", WoE, '\n'), collapse = '')
  sql <- paste0(sql, "     else 0\nend as W", varname, ',\n\n')
  sql
}

toSAS <- function(varname, class, band, WoE)
{
  conditions <- createConditions(varname, class, band, 'SAS')
  SAS <- paste0(paste(c("     if", rep("else if", length(conditions) - 1)), conditions, "then", varname, "=", WoE, ";\n"), collapse = '')
  SAS <- paste0(SAS, "else ", varname, " = 0 ;\n\n")
  SAS
}

save.woe <- function(WoE_result, ...) UseMethod('save.woe')

save.woe.woe.result <- function(WoE_result, outfile)
{
  file_type <- sapply(outfile, function(str) unlist(strsplit(str, '.', fixed = TRUE))[2])


}








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


test <- function(str, ...)
{
  args <- unlist(...)
  return(length(args[[1]]))
}
