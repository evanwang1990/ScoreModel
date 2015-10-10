
#print_trace(deparse(substitute(x)), trace, collapsed_result, method, mode, best_indx, nr, bin_iv)
#自变量名，trace矩阵，分层后的结果，方法，最好的分组，原始的层数,binary_split_iv
print_trace <- function(x_name, trace, collapsed_result, method, mode, best_indx, nr, bin_iv)
{
  if(!is.na(bin_iv)) bin_iv <- round((bin_iv - trace[nr - 1,4]) / bin_iv * 100, 2)
  rownames(trace)[best_indx] <- paste0(rownames(trace)[best_indx],'*')

  trace <- data.frame(trace, row.names = rownames(trace))
  setDT(trace, keep.rownames = T)
  setnames(trace, 'rn', 'Step  ')
  trace[, `:=`(IV             = round(IV, 4),
               IV_decrease     = paste0(round(IV_decrease, 2), '%'),
               X_stat          = round(X_stat, 4),
               C_stat          = round(C_stat, 4),
               Adjust_lift     = round(Adjust_lift, 6),
               Prob.LR_Chi_Sq. = round(Prob.LR_Chi_Sq., 3),
               Prob.z_score.   = round(Prob.z_score., 3),
               Method          = ifelse(Method == 1, 'iv',
                                        ifelse(Method == 2, 'll', 'mo')))]

  cat('Predictor = ', x_name, ', Method = ', method, ', Mode = ', mode, '\n\n', sep = '')
  cat(paste0(rep('-', 100), collapse = ''), '\n')
  print(trace[, c(1, 4, 5, 6, 7, 8, 9, 11, 13, 14), with = F], right = T, row.names = F)
  if(nr > 1 && mode == 'J')
  {
    cat('\nAt last step the iv decrease ', bin_iv, '% from the maximum binary-split-iv.', sep = '')
    if(bin_iv > 5) cat(" The collapse process seems become suboptional.\n") else cat('\n')
  }
  cat('The best collapse is at', trace[[1]][best_indx])
  if(mode == 'J' && trace[best_indx, 9, with = F] > 1e-6) cat(', but it seems not be linear.\n') else cat('.\n')

  if(mode == 'A') setorder(collapsed_result, WoE)
  collapsed_result <- rbind(collapsed_result,
                           data.frame(band       = 'Total',
                                      CntGood    = sum(collapsed_result$CntGood, na.rm = T),
                                      CntBad     = sum(collapsed_result$CntBad, na.rm = T),
                                      CntRec     = sum(collapsed_result$CntRec, na.rm = T),
                                      PctRec     = '100%',
                                      GoodRate   = paste0(round(sum(collapsed_result$CntGood, na.rm = T) / sum(collapsed_result$CntRec, na.rm = T) * 100, 2), '%'),
                                      BadRate    = paste0(round(sum(collapsed_result$CntBad, na.rm = T) / sum(collapsed_result$CntRec, na.rm = T) * 100, 2), '%'),
                                      WoE        = NA,
                                      IV         = sum(collapsed_result$IV),
                                      WoE_barplot= ''
                           ))

  cat("\n\nBest Collapsing Result\n")
  print(collapsed_result)
  cat('\n')
  cat(paste0(rep('=', 100), collapse = ''), '\n\n')
}
