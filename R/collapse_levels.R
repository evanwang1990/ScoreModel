#'
#'library(process)
#'df <- data.frame(x1 = sample(1:1000, 1e4, replace = T),
#'                 x2 = sample(c(letters[1:3], NA), 1e4, replace = T),
#'                 y  = c(sample(c(1, 1, 0), 5e3, replace = T), sample(c(1, 0, 0, 0), 5e3, replace = T)))
#'res <- bestCollapse(c('x1', 'x2'), y, df, 20, method = 'iv', mode = 'J', tracefile = 'trace.Rout', sqlfile = 'sql_code.sql')
#'collapseLevel(x = df$x2, y = df$y, levels = 20, method = 'iv', mode = 'A', minp = 0.05, sourcefile = 'test.R', sqlfile = 'test.sql')

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

print_trace <- function(x_name, trace, method, mode, best_indx, nr, bin_iv)
{
  if(!is.na(bin_iv)) bin_iv <- round((bin_iv - trace[nr - 1,4]) / bin_iv * 100, 1e-3)
  rownames(trace)[best_indx] <- paste0(rownames(trace)[best_indx],'*')
  trace[,c(4, 6, 7, 10, 11, 12)] <- round(trace[,c(4, 6, 7, 10, 11, 12)], 1e-4)
  trace[,5] <- round(trace[,5], 1e-2)
  cat('Predictor = ', x_name, ', Method = ', method, ', Mode = ', mode, '\n', sep = '')
  cat('Mehod: 1->maximum information value; 2->maximum log-likehood; 3->get monotonous.\n\n')
  print(trace)
  if(nr > 1 && mode == 'J')
  {
    cat('\nAt last step the iv decrease ', bin_iv, '% from the maximum binary-split-iv.', sep = '')
    if(bin_iv > 5) cat(" The collapse process seems become suboptional.\n") else cat('\n')
  }
  cat('The best collapse is at', rownames(trace)[best_indx])
  if(mode == 'J' && trace[best_indx, 8] > 1e-6) cat(', but it seems not be linear.\n') else cat('.\n')
  cat("===========================================================\n\n")
}

collapseLevel <- function(x,                                # independent variable
                          y,                                # target
                          levels,                           # initial levels
                          method,                           # collapse methods
                          mode,                             # collapse mode
                          minp = 0.05,                      # min percentage of the sample numbers which one level contains(include missing values)
                          sourcefile,                       # file to store R code
                          sqlfile)
{
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
    x_ <- as.character(x)
  }

  freqMatrix <- as.matrix(table(x_, y))
  if(is.numeric(x))
  {
    labels <- rownames(freqMatrix)
  }else{
    labels <- paste0("'", rownames(freqMatrix), "'")
  }
  rownames(freqMatrix) <- labels

  # deal with the zero cells
  if(any(freqMatrix == 0))
  {
    #collapse levels
    freqMatrix_nonzero <- CollapseZeroCells(freqMatrix, matrix(NA, ncol = 2, nrow = nrow(freqMatrix)), mode = mode)
    # after processing, check if there is still zero cells
    if(any(freqMatrix_nonzero$freqMatrix == 0))
    {
      warning("There are zeros in some cells!\n")
      return(list("need_mo" = FALSE,
                  "main"    = data.frame(var         = deparse(substitute(x)),
                                         class       = class(x),
                                         NAs         = sum(is.na(x))/length(x),
                                         iv          = NA,
                                         levels      = sum(groups != " "),
                                         linear      = ifelse(mode == 'J', trace[best_indx, 7] < 1e-6, NA),
                                         suboptional = NA,
                                         method      = method,
                                         detail      = 'zero-cell')))
    }

    #update labels of freqMatrix
    new_labels <- GetGroups(labels, freqMatrix_nonzero$label.trace[!is.na(freqMatrix_nonzero$label.trace[,1]),1], freqMatrix_nonzero$label.trace[!is.na(freqMatrix_nonzero$label.trace[,2]),2])
    freqMatrix <- freqMatrix_nonzero$freqMatrix
    rownames(freqMatrix) <- new_labels[new_labels != " "]
  }

  #check nrow
  nr <- nrow(freqMatrix)
  if(nr > 1)
  {
    trace <- matrix(nrow = nr,
                    ncol = 13,
                    dimnames = list(paste("Step", 0:(nr - 1)),
                                    c('Left', 'Right', 'minCount', 'IV', 'IV decrease %', 'X_stat', 'C_stat', 'Adjust lift', 'Log likehood', 'Prob(x > LR_Chi_Sq)', 'Z_score of log odds ratio', 'Prob(z_score = 0)', 'Method')))
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 10] <- 1 - pchisq(-2 * (trace[-1, 9] - trace[1, 9]), 1:(nr - 1))
    trace[-1, 12] <- 1 - 2 * abs(pnorm(trace[-1, 11]) - 0.5)
    bin_iv <- trace[1, 12]
    trace[1, 12] <- 1
    best_indx <- max(min(which.max(trace[, 10] < 0.05) - 1, which.max(trace[, 12] < 0.05) - 1), which.max(trace[, 3] >= minp * length(x)))
    print_trace(deparse(substitute(x)), trace, method, mode, best_indx, nr, bin_iv)

    #string code
    if(best_indx == 1)
      groups <- rownames(freqMatrix)
    else
      groups <- GetGroups(rownames(freqMatrix), trace[2:best_indx, 1], trace[2:best_indx, 2])
    stringcode <- stringCode(x, y, x_, deparse(substitute(x)), groups, !missing(sqlfile), method)

    #output
    writeLines(stringcode[1], sourcefile, sep = ',\n')
    if(!missing(sqlfile)) writeLines(stringcode[2], sqlfile)


    return(list("need_mo" = (mode == 'J' && method != 'mo' && trace[best_indx, 8] > 1e-6),
                "main"    = data.frame(var         = deparse(substitute(x)),
                                       class       = class(x),
                                       NAs         = sum(is.na(x))/length(x),
                                       iv          = round(trace[best_indx, 4], 1e-4),
                                       levels      = sum(groups != " "),
                                       linear      = ifelse(mode == 'J', trace[best_indx, 8] < 1e-6, NA),
                                       suboptional = ifelse(mode == 'J', round((bin_iv - trace[nr - 1,4]) / bin_iv * 100, 1e-3), NA),
                                       method      = method,
                                       detail      = 'normal',
                                       row.names   = NULL)))

  }else if(any(is.na(x))){  #the one level variables make sense only when they have missing values(just as two levels)
    stringcode <- stringCode(x, y, x_, deparse(substitute(x)), paste0("'",rownames(freqMatrix),"'"), !missing(sqlfile), method)
    writeLines(stringcode[1], sourcefile, sep = ',\n')
    if(!missing(sqlfile)) writeLines(stringcode[2], sqlfile)
    return(list("need_mo" = FALSE,
                "main"    = data.frame(var         = deparse(substitute(x)),
                                       class       = class(x),
                                       NAs         = sum(is.na(x))/length(x),
                                       iv          = 0,
                                       levels      = sum(groups != " "),
                                       linear      = ifelse(mode == 'J', trace[best_indx, 8] < 1e-6, NA),
                                       suboptional = 0,
                                       method      = method,
                                       detail      = 'one-level')))
  }
}

bestCollapse <- function(vars,
                         target,
                         dataset,
                         max.levels,
                         method = c('iv', 'll', 'mo'),
                         mode = 'J',
                         tracefile,
                         sourcefile,
                         sqlfile)
{
  attach(dataset, warn.conflicts = F)
  # check the parameters
  if(any(is.na(target)))                                        #target
    stop("There are NAs in target!\n")
  out_vars <- vars[!vars %in% names(dataset)]                   #vars
  if(length(out_vars) != 0)
    stop(paste(out_vars, sep = ", "), " are not in ", deparse(substitute(dataset)), " !\n")
  method <- match.arg(method)                                   #method
  if(!missing(tracefile))                                       #tracefile
  {
    tracefile <- file(tracefile, open = 'wt')
    sink(file = tracefile, append = T, type = 'output')
  }
  if(missing(sourcefile))                                       #sourcefile
  {
    warning("`sourcefile` is missing, default file named 'woe_code.R' is set.\n")
    sourcefile <- 'woe_code.R'
  }
  sourcefile <- file(sourcefile, open = 'wt')
  writeLines("library(data.table)", sourcefile)
  writeLines(paste0("setDT(", deparse(substitute(dataset)), ")"), sourcefile)
  writeLines(paste0(deparse(substitute(dataset)), "[ , `:=`("), sourcefile)
  if(!missing(sqlfile))                                         #sqlfile
    sqlfile <- file(sqlfile, open = 'wt')

  for(var in vars)
  {
    res <- do.call(collapseLevel, list(as.name(var), target, max.levels, method, mode, sourcefile, sqlfile))
    if(!exists('main'))
    {
      main <- res$main
    }else{
      main <- rbind(main, res$main)
    }
    if(res$need_mo)
    {
      res <- do.call(collapseLevel, list(as.name(var), target, max.levels, 'mo', mode, sourcefile, sqlfile))
      main <- rbind(main, res$main)
    }
  }

  writeLines(")]", sourcefile, sep = ';\n')
  close(sourcefile)
  if(!missing(tracefile)) close(tracefile)
  if(!missing(sqlfile)) close(sqlfile)
  main
}
