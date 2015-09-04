#include <Rcpp.h>
#include "calculate_stat.h"
using    namespace Rcpp;
using    namespace std;

double binary_split(NumericMatrix freqMatrix);
//test:
//tmp <- sample(20:100, 10, replace = T)
//tmp <- matrix(c(tmp, round(tmp * (seq(1, 3, length.out = 10) + 0.5 * runif(1:10)) + sample(1:10, 10))), ncol = 2)
//rownames(tmp) <- letters[1:10]
//plot(tmp[,2]/tmp[,1])


//[[Rcpp::export]]
NumericMatrix Collapse(NumericMatrix freqMatrix, NumericMatrix trace, int row_indx = 1, String method = "iv", String mode = "J")
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad  = wrap(freqMatrix.column(1));
  size_t n           = good.size();
  double min_delta   = R_PosInf, delta_;
  int left, right;

  //calculate initial log-likehood and max-binary-informance-value
  if(row_indx == 1)
  {
    trace(0, 8)  = cal_ll(good, bad);
    trace(0, 10) = binary_split(freqMatrix);
    trace(0, 3) = cal_iv(good / sum(good), bad / sum(bad));
  }

  //find then best collapse
  if(mode == "J")
  {
    for(size_t i = 0; i < n - 1; ++i)
    {
      delta_ = delta(good, bad, i, i + 1, method);
      //when using method 'mo',if it's already monotonous
      //or there's no levels to collapse use method 'iv'
      if(delta_ < min_delta)
      {
        min_delta = delta_;
        left      = i;
        right     = i + 1;
      }
    }
  }
  else
  {
    for(size_t i = 0; i < n - 1; ++i)
    {
      for(size_t j = i + 1; j < n; ++j)
      {
        delta_ = delta(good, bad, i, j, method);
        if(delta_ < min_delta)
        {
          min_delta = delta_;
          left      = i;
          right     = j;
        }
      }
    }
  }

  if(min_delta == 9999.0) return(Collapse(freqMatrix, trace, row_indx, "iv", "J"));

  //get new freqmatrix
  NumericMatrix freqMatrix_ = combine(freqMatrix, left, right);
  NumericVector new_good    = wrap(freqMatrix_.column(0));
  NumericVector new_bad     = wrap(freqMatrix_.column(1));

  //trace
  trace(row_indx, 0) = left;
  trace(row_indx, 1) = right;
  trace(row_indx, 9) = cal_log_odds_ratio_zscore(good, bad, left, right);
  trace(row_indx, 7) = cal_ll(new_good, new_bad);
  new_good           = new_good / sum(new_good);
  new_bad            = new_bad / sum(new_bad);
  trace(row_indx, 2) = cal_iv(new_good, new_bad);
  trace(row_indx, 3) = (trace(row_indx, 2) - trace(0, 3)) * 100 / trace(0, 3);
  trace(row_indx, 4) = cal_x_stat(new_good, new_bad);
  trace(row_indx, 5) = mode == "J" ? cal_c_stat(new_good, new_bad):NA_REAL;
  trace(row_indx, 6) = mode == "J" ? (trace(row_indx, 4) - trace(row_indx, 5)) / (trace(row_indx, 5) * (n - 2)):NA_REAL;
  trace(row_indx, 11) = method == "iv" ? 1:(method == "ll" ? 2:3);

  if(freqMatrix_.nrow() == 2) return(trace);
  return(Collapse(freqMatrix_, trace, ++ row_indx, method, mode));
}

//[[Rcpp::export]]
StringVector GetGroups(StringVector labels, NumericVector left, NumericVector right)
{
  StringVector labels_ = clone(labels);
  for(int i = 0; i < left.size(); ++i)
  {
    combineLabels(labels_, left[i], right[i]);
  }
  return(labels_);
}

double binary_split(NumericMatrix freqMatrix)
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad  = wrap(freqMatrix.column(1));
  size_t nr          = freqMatrix.nrow();
  double max_iv      = R_NegInf, iv;

  for(size_t i = 1; i < nr - 1; ++i)
  {
    double good1 = 0, good2 = 0, bad1 = 0, bad2 = 0;
    for(size_t j = 0; j < nr; ++j)
    {
      if(j < i)
      {
        good1 += good[j];
        bad1  += bad[j];
      }
      else
      {
        good2 += good[j];
        bad2  += bad[j];
      }
    }

    iv = (good1/(good1 + good2) - bad1/(bad1 + bad2)) * (log(good1/(good1 + good2)) - log(bad1/(bad1 + bad2))) +
      (good2/(good1 + good2) - bad2/(bad1 + bad2)) * (log(good2/(good1 + good2)) - log(bad2/(bad1 + bad2)));
    max_iv = max(max_iv, iv);
  }
  return(max_iv);
}
