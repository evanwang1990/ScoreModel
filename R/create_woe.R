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
