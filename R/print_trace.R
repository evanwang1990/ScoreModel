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
