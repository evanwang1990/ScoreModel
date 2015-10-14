catLog <- function(WoE_result) UseMethod('catLog')

catLog.woe.result <- function(WoE_result)
{
  cat('Predictor = ', WoE_result$summary$var, ', Method = ', WoE_result$summary$method, ', Mode = ', WoE_result$summary$mode, '\n\n', sep = '')
  cat(paste0(rep('-', 100), collapse = ''), '\n\n')

  # print trace
  trace <- WoE_result$trace
  if(!is.null(trace))
  {
    trace <- data.frame(trace, row.names = rownames(trace), stringsAsFactors = F)
    setDT(trace, keep.rownames = T)
    setnames(trace, 'rn', paste0(c('Step', rep(' ', max(nchar(trace[[1]])) - 4)), collapse = ''))
    trace[, `:=`(IV              = round(IV, 4),
                 IV_decrease     = paste0(round(IV_decrease, 2), '%'),
                 X_stat          = round(X_stat, 4),
                 C_stat          = round(C_stat, 4),
                 Adjust_lift     = round(Adjust_lift, 6),
                 Prob.LR_Chi_Sq. = round(Prob.LR_Chi_Sq., 3),
                 Prob.z_score.   = round(Prob.z_score., 3),
                 Method          = ifelse(Method == 1, 'max_iv',
                                          ifelse(Method == 2, 'max_likehood', 'linear')))]
    cat('Trace of collapsing levels:\n\n')
    print(trace[, c(1, 4, 5, 6, 7, 8, 9, 11, 13, 14), with = F], right = T, row.names = F)
    if (WoE_result$summary$levels > 1 && WoE_result$summary$mode == "J")
    {
      cat('\nAt last step the iv decrease ', binary_IV, '% from the maximum binary-split-iv.', sep = '')
      if (binary_IV > 5) cat(" The collapse process seems become suboptional.\n") else cat('\n')
    }
    cat(paste0(rep('-', 100), collapse = ''), '\n\n')
  }

  #print WoE_result
  cat('The result of WoE:\n')
  print(WoE_result$detail, right = T, row.names = F)
  cat('\n')
  if(!WoE_result$summary$is.linear) cat('## It seems NOT to be linear.\n')
  else cat('## It is linear.\n')
  cat(paste0(rep('=', 100), collapse = ''), '\n\n\n\n')
}

catLog.woe.results <- function(WoE_result)
{
  lapply(WoE_result, catLog.woe.result)
}
