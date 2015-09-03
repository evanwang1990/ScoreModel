#include <Rcpp.h>
#include "calculate_stat.h"
using namespace Rcpp;
using namespace std;


//[[Rcpp::export]]
List sub_collapse(NumericMatrix freqMatrix, String method = "iv", String mode = "J")
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad  = wrap(freqMatrix.column(1));
  size_t n = good.size();
  double min_delta = R_PosInf, delta_;
  int left, right;

  //find then best collapse
  if(mode == "J")
  {
    for(size_t i = 0; i < n - 1; ++i)
    {
      delta_ = delta(good, bad, i, i + 1, method);
      //when using method 'mo',if it's already monotonous
      //or there's no levels to collapse use method 'iv'
      if(delta_ == 9999.0) return(sub_collapse(freqMatrix));
      if(delta_ < min_delta)
      {
        min_delta = delta_;
        left = i; right = i + 1;
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
          left = i; right = j;
        }
      }
    }
  }

  //get new freqmatrix
  NumericMatrix freqMatrix_ = combine(freqMatrix, left, right);
  NumericVector new_good    = wrap(freqMatrix_.column(0));
  NumericVector new_bad     = wrap(freqMatrix_.column(1));

  //get statistics in each step
  double lo_zscore = cal_log_odds_ratio_zscore(good, bad, left, right);
  double ll        = cal_ll(new_good, new_bad);
  new_good         = new_good / sum(new_good);
  new_bad          = new_bad / sum(new_bad);
  double iv        = cal_iv(new_good, new_bad);
  double x_stat    = cal_x_stat(new_good, new_bad);
  double c_stat    = method == "J" ? cal_c_stat(new_good, new_bad):NA_REAL;
  double adj_lift  = method == "J" ? (x_stat - c_stat) / (c_stat * (n - 2)):NA_REAL;

  return(List::create(Named("freqMatrix") = freqMatrix_,
                      Named("stat") = NumericVector::create(left, right, iv, x_stat, c_stat, adj_lift, ll, NA_REAL, lo_zscore, NA_REAL, NA_REAL)));
}

//[[Rcpp::export]]
double binary_split(NumericMatrix freqMatrix)
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad  = wrap(freqMatrix.column(1));
  size_t nr          = freqMatrix.nrow();
  double max_iv = R_NegInf, iv;
  for(size_t i = 1; i < nr - 1; ++i)
  {
    double good1 = 0; double good2 = 0; double bad1 = 0; double bad2 = 0;
    for(size_t j = 0; j < nr; ++j)
    {
      if(j < i)
      {
        good1 += good[j];
        bad1 += bad[j];
      }
      else
      {
        good2 += good[j];
        bad2 += bad[j];
      }
    }

    iv = (good1/(good1 + good2) - bad1/(bad1 + bad2)) * (log(good1/(good1 + good2)) - log(bad1/(bad1 + bad2))) +
         (good2/(good1 + good2) - bad2/(bad1 + bad2)) * (log(good2/(good1 + good2)) - log(bad2/(bad1 + bad2)));
    max_iv = max(max_iv, iv);
    cout<<"i = "<<i<<" iv"<<iv<<endl;
  }
  return(max_iv);
}

//[[Rcpp::export]]
List collapse(NumericMatrix freqMatrix, String method = "iv", String mode = "J")
{
  size_t nr = freqMatrix.nrow();
  List dimnames = freqMatrix.attr("dimnames");
  StringVector labels = wrap(dimnames[0]);
  labels = clone(labels);

  //the columns are:left, right, iv, x_stat, c_stat, adjust-lift, log likehood, likehood ratio chisq, log_odds ratio z_score, max binary split IV
  NumericMatrix trace(nr - 1, 10);
  trace.fill(NA_REAL);
  //store labels while collapsing levels
  StringMatrix traceLabels(nr - 1, nr);

  //initialization
  NumericVector good0 = wrap(freqMatrix.column(0));
  NumericVector bad0  = wrap(freqMatrix.column(1));
  double ll0 = cal_ll(good0, bad0);
  good0 = good0 / sum(good0);
  bad0  = bad0 / sum(bad0);
  double iv0 = cal_iv(good0, bad0);
  double x_stat0 = cal_x_stat(good0, bad0);
  double c_stat0 = cal_c_stat(good0, bad0);
  trace(0,2) = iv0;
  trace(0,3) = x_stat0;
  trace(0,4) = c_stat0;
  trace(0,6) = ll0;

  traceLabels.row(0) = labels;

  //collapse
  List sub_collapse_result;
  sub_collapse_result = sub_collapse(freqMatrix, method, mode);
  for(size_t i = 1; i < nr - 1; ++i)
  {
    NumericMatrix freqMatrix_ = sub_collapse_result["freqMatrix"];
    NumericVector stat = sub_collapse_result["stat"];
    trace.row(i) = stat;
    //NumericVector loc = sub_collapse_result["loc"];
    //trace(i,0)  = loc[0];
   // trace(i,1)  = loc[1];
   // trace(i,2)  = sub_collapse_result["iv"];
   // trace(i,3)  = sub_collapse_result["x_stat"];
    //monotony(c_stat) makes no sense in "All pairs" method
   // if(mode == "J")
  //  {
  //    trace(i,4)  = sub_collapse_result["c_stat"];
  //    trace(i,5)  = (trace(i,3) - trace(i,4)) / (trace(i, 4) * (nr - i - 1));
  //  }
  //  trace(i,6)  = sub_collapse_result["ll"];
  //  trace(i,7)  = -2 * (trace(i,6) - ll0);
  //  trace(i,8)  = sub_collapse_result["lo_zscore"];
  //  if(i == nr - 2)
  //  {
 //     //calculate max binary split iv
  //    double max_binary_iv = binary_split(freqMatrix);
  //    trace(i,9) = (max_binary_iv - trace(i,2)) / max_binary_iv * 100;
  //  }

    combineLabels(labels, stat[0], stat[1]);
    traceLabels.row(i) = labels;

    sub_collapse_result = sub_collapse(freqMatrix_, method, mode);
  }

  return(List::create(Named("trace") = trace,
                      Named("traceLabels") = traceLabels));
}

//test:
//tmp <- sample(20:100, 10, replace = T)
//tmp <- matrix(c(tmp, round(tmp * (seq(1, 3, length.out = 10) + 0.5 * runif(1:10)) + sample(1:10, 10))), ncol = 2)
//rownames(tmp) <- letters[1:10]
//plot(tmp[,2]/tmp[,1])
