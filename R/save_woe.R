#注意：变量长度太长时，需要截断加密，制作对应表map

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
  sql <- paste0(sql, "     else 0\nend as W", varname, ',\n')
  sql
}

toSAS <- function(varname, class, band, WoE)
{
  conditions <- createConditions(varname, class, band, 'SAS')
  SAS <- paste0(paste0(c("     if ", rep("else if ", length(conditions) - 1)), conditions, " then W", varname, " = ", WoE, " ;\n"), collapse = '')
  SAS <- paste0(SAS, "else W", varname, " = 0 ;\n")
  SAS
}

save.woe <- function(WoE_results, outfile)
{
  if (!is.list(outfile))
  {
    if (is.character(outfile)) outfile <- as.list(outfile)
    else stop("outfile should be a list or a character vector")
  }
  file_type <- sapply(outfile, function(str) unlist(strsplit(str, '.', fixed = TRUE))[2])
  if (any(!file_type %in% c('SAS', 'sql'))) stop("the format of the file to save in not supported")
  sapply(outfile, function(x) if (file.exists(x)) file.remove(x))

  for (WoE_result in WoE_results)
  {
    for (i in 1:length(outfile))
    {
      stringcode <- ifelse(file_type[i] == 'sql',
                           toSql(WoE_result$summary$var, WoE_result$summary$class, WoE_result$detail$band[1:(nrow(WoE_result$detail) - 1)], WoE_result$detail$WoE[1:(nrow(WoE_result$detail) - 1)]),
                           toSAS(WoE_result$summary$var, WoE_result$summary$class, WoE_result$detail$band[1:(nrow(WoE_result$detail) - 1)], WoE_result$detail$WoE[1:(nrow(WoE_result$detail) - 1)]))
      write(stringcode, outfile[[i]], append = T, sep = '\n')
    }
  }
}
