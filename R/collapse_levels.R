#'
#'library(process)
#'df <- data.frame(x1 = sample(1:1000, 1e4, replace = T),
#'                 x2 = sample(c(letters[1:3], NA), 1e4, replace = T),
#'                 y  = c(sample(c(1, 1, 0), 5e3, replace = T), sample(c(1, 0, 0, 0), 5e3, replace = T)))
#'res <- bestCollapse(c('x1', 'x2'), y, df, 20, method = 'iv', mode = 'J', tracefile = 'trace.Rout', sqlfile = 'sql_code.sql')
#'collapseLevel(x = df$x2, y = df$y, levels = 20, method = 'iv', mode = 'A', minp = 0.05, sourcefile = 'test.R', sqlfile = 'test.sql')
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
  }else{
    mode <- 'J'
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
                                         levels      = NA,
                                         linear      = NA,
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
                                    c('Left', 'Right', 'minCount', 'IV', 'IV_decrease', 'X_stat', 'C_stat', 'Adjust_lift', 'Log_likehood', 'Prob(LR_Chi_Sq)', 'Z_score_of_log_odds_ratio', 'Prob(z_score)', 'Method')))
    trace <- Collapse(freqMatrix, trace, 1, method, mode)

    #choose the best collapse
    trace[-1, 10] <- 1 - pchisq(-2 * (trace[-1, 9] - trace[1, 9]), 1:(nr - 1))
    trace[-1, 12] <- 1 - 2 * abs(pnorm(trace[-1, 11]) - 0.5)
    bin_iv <- trace[1, 12]
    trace[1, 12] <- 1
    best_indx <- max(min(which.max(trace[, 10] < 0.05) - 1, which.max(trace[, 12] < 0.05) - 1), which.max(trace[, 3] >= minp * length(x)))

    #get collapsed result----
    if(best_indx == 1){
      collapsed_result <- freqMatrix
    }else{
      collapsed_result <- combineResults(freqMatrix, trace[2:best_indx, 1], trace[2:best_indx, 2])
    }
    collapsed_result <- data.table(collapsed_result, keep.rownames = T)
    setnames(collapsed_result, c('band', 'CntGood', 'CntBad'))
    collapsed_result[, `:=`(band         = band(x, x_, band),
                            CntRec       = CntGood + CntBad)
                     ][, `:=`(PctRec     = paste0(round(CntRec / length(x) * 100, 2), '%'),
                              GoodRate   = paste0(round(CntGood / CntRec * 100, 2), '%'),
                              BadRate    = paste0(round(CntBad / CntRec * 100, 2), '%'),
                              WoE        = round(log(CntGood / sum(y == 0)) - log(CntBad / sum(y == 1)), 4),
                              IV         = round((CntGood / sum(y == 0) - CntBad / sum(y == 1)) * (log(CntGood / sum(y == 0)) - log(CntBad / sum(y == 1))), 4))]
    if(any(is.na(x))) #add missing values
    {
      NA_good <- sum(is.na(x) & y == 0)
      NA_bad  <- sum(is.na(x) & y == 1)
      NA_tot <- NA_good + NA_bad

      collapsed_result <- rbind(collapsed_result,
                               data.frame(band       = 'missing',
                                           CntGood    = NA_good,
                                           CntBad     = NA_bad,
                                           CntRec     = NA_tot,
                                           PctRec     = paste0(round(NA_tot / length(x) * 100, 2), '%'),
                                           GoodRate   = paste0(round(NA_good / NA_tot * 100, 2), '%'),
                                           BadRate    = paste0(round(NA_bad / NA_tot * 100, 2), '%'),
                                           WoE        = ifelse(NA_good * NA_bad == 0, NA, round(log(NA_good / sum(y == 0)) - log(NA_bad / sum(y == 1)), 4)),
                                           IV         = ifelse(NA_good * NA_bad == 0, NA, round((NA_good / sum(y == 0) - NA_bad / sum(y == 1)) * (log(NA_good / sum(y == 0)) - log(NA_bad / sum(y == 0))), 4))
                                           ))
    }

    collapsed_result[, WoE_barplot := barplot.woe(WoE)]

    #print trace----
    print_trace(deparse(substitute(x)), trace, collapsed_result, method, mode, best_indx, nr, bin_iv)


    #stringcode <- stringCode(x, y, x_, deparse(substitute(x)), groups, !missing(sqlfile), method)

    #output
    #writeLines(stringcode[1], sourcefile, sep = ',\n')
    #if(!missing(sqlfile)) writeLines(stringcode[2], sqlfile)


    return(list("need_mo" = (mode == 'J' && method != 'mo' && trace[best_indx, 8] > 1e-6),
                "main"    = data.frame(var         = deparse(substitute(x)),
                                       class       = class(x),
                                       NAs         = sum(is.na(x))/length(x),
                                       iv          = round(trace[best_indx, 4], 1e-4),
                                       levels      = nrow(collapsed_result) - any(is.na(x)),
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
                                       levels      = 1,
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
