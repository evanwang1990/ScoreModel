stringCode <- function(x, y, x_, x_name, groups, SQLcode)
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
  code <- paste0(x_name, "_w = ", paste0("recode(", paste(x_name, code, "as.factor.result = F", sep = ','), ")"))

  #string SQL code
  if(SQLcode)
  {
    sql_code <- paste('\twhen', range_sql, 'then', woe[!is.na(names(woe))], collapse = "\n")
    if(any(is.na(x))) sql_code <- paste(sql_code, '\n\telse', tail(woe,1))
    sql_code <- paste('case', sql_code, paste('as', paste0(x_name, '_w'), ';'), sep = '\n')
    return(c(code, sql_code))
  }
  return(code)
}

print_trace <- function(x_name, trace, method, mode, best_indx)
{
  bin_iv <- round((trace[1,11] - tail(trace[,3],1)) / trace[1,11] * 100, 1e-3)
  rownames(trace)[best_indx] <- paste0(rownames(trace)[best_indx],'*')
  trace[,c(3, 5, 6, 9, 10, 11)] <- round(trace[,c(3, 5, 6, 9, 10, 11)], 1e-4)
  trace[,4] <- round(trace[,4], 1e-2)
  cat('Predictor = ', x_name, ', Method = ', method, ', Mode = ', mode, '\n', sep = '')
  cat('Mehod: 1->maximum information value;2->maximum log-likehood;3->get monotonous\n\n')
  print(trace[-1,])
  cat('\nAt last step the iv decrease ', bin_iv, '% from the maximum binary-split-iv.', sep = '')
  if(bin_iv > 5) cat(" The collapse process seems become suboptional.\n") else cat('\n')
  cat('The best collapse is at', rownames(trace)[best_indx])
  if(mode == 'J' && trace[best_indx, 7] > 1e-6) cat(', but it seems not be linear.\n') else cat('\n')
  cat("===========================================================\n\n")
}

#x <- sample(c(1:1000, rep(NA, 100)), 1e4, replace = T)
#x <- sample(c(letters[1:10], NA), 1e4, replace = T)
#x <- sample(c(letters[1:10]), 1e4, replace = T)
#y <- sample(c(1, 0), 1e4, replace = T)
collapseLevel <- function(x,                                # independent variable
                          y,                                # target
                          levels,                           # initial levels
                          method = c('iv', 'll', 'mo'),     # collapse methods
                          mode = 'J',                       # collapse mode
                          do.trace = T,                     # need output collapsing step?
                          sourcefile,                       # file to store R code
                          sqlfile)
{
  if(missing(sourcefile)) stop("`sourcefile` is missing!\n")
  # deal with NAs
  if(any(is.na(y))) stop("There are NAs in target variable!\n")

  # deal with `method` and `mode`
  method <- match.arg(method)
  if(is.character(x) || (is.factor(x) && !is.ordered(x)))
  {
    if(method == 'mo')
    {
      warning("Independent variable is character vector or unordered vector, 'mo' method is not appreciate\n The method is modified into 'iv' automatically.")
      method <- 'iv'
    }
    mode <- 'A'
  }

  if(is.numeric(x))
  {
    x_ <- cut(rank(x, ties.method = 'min', na.last = "keep"), breaks = levels, labels = F)
  }else{
    x_ <- x
  }


  freqMatrix <- as.matrix(table(x_, y))
  # check if there are zeros in cells
  if(any(freqMatrix == 0)) stop("There are zero cells!\n")


  #check nrow
  nr <- nrow(freqMatrix)
  if(nrow(freqMatrix) > 2)
  {
    trace <- matrix(nrow = nr - 1,
                    ncol = 12,
                    dimnames = list(paste("Step", 0:(nr - 2)),
                                    c('Left', 'Right', 'IV', 'IV decrease %', 'X_stat', 'c_stat', 'Adjust lift', 'Log likehood', 'Prob(x > LR_Chi_Sq)', 'Z_score of log odds ratio', 'Prob(z_score = 0)', 'Method')))
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 9] <- 1 - pchisq(-2 * (trace[-1, 8] - trace[1, 9]), (nr - 2) : 1)
    trace[-1, 11] <- 1 - 2 * abs(pnorm(trace[-1, 10]) - 0.5)
    best_indx <- sum((trace[-1, 9] < 0.1) + (trace[-1, 11] < 0.1)  == 0) + 1
    print_trace(deparse(substitute(x)), trace, method, mode, best_indx)

    #string code
    if(is.numeric(x))
    {
      labels <- rownames(freqMatrix)
    }else{
      labels <- paste0("'", rownames(freqMatrix), "'")
    }
    groups <- GetGroups(labels, trace[2:best_indx, 1], trace[2:best_indx, 2])
    SQLcode <- !missing(sqlfile)
    stringcode <- stringCode(x, y, x_, deparse(substitute(x)), groups, SQLcode)

    #output
    writeLines(stringcode[1], sourcefile)
    if(SQLcode) writeLines(stringcode[2], sqlfile)


    return(list("need_mo" = (mode == 'J' && method != 'mo' && trace[best_indx, 7] > 1e-6),
                "main"    = data.frame(var         = deparse(substitute(x)),
                                       iv          = round(trace[best_indx, 3], 1e-4),
                                       method      = method,
                                       suboptional = round((trace[1,11] - tail(trace[,3],1)) / trace[1,11] * 100, 1e-3),
                                       row.names   = NULL)))

  }
}
