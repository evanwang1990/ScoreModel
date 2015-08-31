#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
inline double cal_c_stat(NumericVector good, NumericVector bad)
{
  double res = 0;
  size_t n = good.size();
  for(size_t i = 0; i < n; ++i)
  {
    for(size_t j = i; j < n; ++j)
    {
      res += (i == j ? 0.5:1) * good[i] * bad[j];
    }
  }

  res = max(res, 1-res);
  return(res);
}

//[[Rcpp::export]]
inline double cal_x_stat(NumericVector good, NumericVector bad)
{
  double x_stat = 0;
  size_t n = good.size();
  for(size_t i = 0; i < n; ++i)
  {
    for(size_t j = i + 1; j < n; ++j)
    {
      x_stat += abs(good[i] * bad[j] - good[j] * bad[i]);
    }
  }

  x_stat = (1 + x_stat) / 2;
  return(x_stat);
}

//[[Rcpp::export]]
double cal_iv(NumericVector good, NumericVector bad)
{
  double iv = 0;
  size_t n = good.size();
  for(size_t i = 0; i < n; ++i)
  {
    iv += (good[i] - bad[i]) * log(good[i] / bad[i]);
  }
  return(iv);
}

//[[Rcpp::export]]
double cal_ll(NumericVector good, NumericVector bad)
{
  double ll = 0;
  int tot;
  size_t n = good.size();

  for(size_t i = 0; i < n; ++i)
  {
    tot = good[i] + bad[i];
    ll += good[i] * log(good[i] / tot) + bad[i] * log(bad[i] / tot);
  }

  return(ll);
}

//[[Rcpp::export]]
inline double cal_log_odds_ratio_zscore(NumericVector good, NumericVector bad, int i, int j)
{
  double z_score;
  z_score = log((good[i] * bad[j]) / (good[j] * bad[i])) / sqrt(1 / good[i] + 1 / good[j] + 1 / bad[i] + 1 / bad[j]);
  return(z_score);
}

//[[Rcpp::export]]
inline NumericMatrix combine(NumericMatrix freqMatrix, int i, int j)
{
  size_t nr = freqMatrix.nrow();
  size_t nc = freqMatrix.ncol();
  NumericMatrix freqMatrix_ = NumericMatrix(nr - 1, nc);

  for(size_t k = 0; k < nr - 1; ++k)
  {
    if(k == i)
    {
      freqMatrix_.row(k) = freqMatrix.row(i) + freqMatrix.row(j);
    }else if(k < j)
    {
      freqMatrix_.row(k) = freqMatrix.row(k);
    }else if(k >= j)
    {
      freqMatrix_.row(k) = freqMatrix.row(k + 1);
    }
  }

  return(freqMatrix_);
}

//[[Rcpp::export]]
double delta(NumericVector good, NumericVector bad, int i, int j, String method = "iv")
{
  double delta_;
  if(method == "iv")
  {
    double sum_good = sum(good);
    double sum_bad = sum(bad);
    delta_ = (good[i] / sum_good - bad[i] / sum_bad) * (log(good[i] / sum_good) - log(bad[i] / sum_bad)) +
            (good[j] / sum_good - bad[j] / sum_bad) * (log(good[j] / sum_good) - log(bad[j] / sum_bad)) -
            ((good[i] + good[j]) / sum_good - (bad[i] + bad[j]) / sum_bad) * (log((good[i] + good[j]) / sum_good) - log((bad[i] + bad[j]) / sum_bad));
  }
  else if(method == "ll")
  {
    double new_good = good[i] + good[j];
    double new_bad = bad[i] + bad[j];
    delta_ = good[i] * log(good[i] / (good[i] + bad[i])) + bad[i] * log(bad[i] / (good[i] + bad[i])) +
            good[j] * log(good[j] / (good[j] + bad[j])) + bad[j] * log(bad[j] / (good[j] + bad[j])) -
            (new_good * log(new_good / (new_good + new_bad)) + new_bad * log(new_bad / (new_good + new_bad)));
  }
  // make sure that only numeric x can use "mo" method!!
  else if(method == "mo")
  {
    int n = good.size() - 1;
    NumericVector new_good(n);
    NumericVector new_bad(n);
    for(int k = 0; k < n; ++k)
    {
      if(k == i)
      {
        new_good[k] = good[i] + bad[j];
        new_bad[k]  = bad[i] + bad[j];
      }
      else if(k <= j)
      {
        new_good[k] = good[k];
        new_bad[k]  = bad[k];
      }
      else
      {
        new_good[k] = good[k + 1];
        new_bad[k] = bad[k + 1];
      }
    }

    new_good = new_good / sum(new_good);
    new_bad = new_bad / sum(new_bad);

    delta_ = cal_x_stat(new_good, new_bad) - cal_c_stat(new_good, new_bad);
  }

  return(delta_);
}

//[[Rcpp::export]]
List sub_collapse(NumericMatrix freqMatrix, String method = "iv", String mode = "J")
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad = wrap(freqMatrix.column(1));
  size_t n = good.size();
  double min_delta = R_PosInf;
  double delta_;
  IntegerVector loc(2);

  //find then best collapse
  if(mode == "J")
  {
    for(size_t i = 0; i < n - 1; ++i)
    {
      delta_ = delta(good, bad, i, i + 1, method);
      if(delta_ < min_delta)
      {
        min_delta = delta_;
        loc[0] = i;
        loc[1] = i + 1;
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
          loc[0] = i;
          loc[1] = j;
        }
      }
    }
  }


  //get new freqmatrix
  NumericMatrix freqMatrix_ = combine(freqMatrix, loc[0], loc[1]);
  NumericVector new_good = wrap(freqMatrix_.column(0));
  NumericVector new_bad = wrap(freqMatrix_.column(1));

  //get indexes in each step
  double lo_zscore = cal_log_odds_ratio_zscore(good, bad, loc[0], loc[1]);
  double ll = cal_ll(new_good, new_bad);
  new_good = new_good / sum(new_good);
  new_bad = new_bad / sum(new_bad);
  double iv = cal_iv(new_good, new_bad);
  double x_stat = cal_x_stat(new_good, new_bad);
  double c_stat = cal_c_stat(new_good, new_bad);


  return(List::create(Named("freqMatrix") = freqMatrix_,
                      Named("loc") = loc,
                      Named("iv") = iv,
                      Named("x_stat") = x_stat,
                      Named("c_stat") = c_stat,
                      Named("ll") = ll,
                      Named("lo_zscore") = lo_zscore));
}

//[[Rcpp::export]]
NumericMatrix collapse(NumericMatrix freqMatrix, String method = "iv", String mode = "J")
{
  size_t nr = freqMatrix.nrow();
  //the columns are:left, right, iv, x_stat, c_stat, adjust-lift, log likehood, likehood ratio chisq, log_odds ratio z_score
  NumericMatrix trace(nr - 1, 9);
  trace.fill(NA_REAL);

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

  //collapse
  List sub_collapse_result;
  sub_collapse_result = sub_collapse(freqMatrix, method, mode);
  for(size_t i = 1; i < nr - 1; ++i)
  {
    NumericMatrix freqMatrix_ = sub_collapse_result["freqMatrix"];
    NumericVector loc = sub_collapse_result["loc"];
    trace(i,0)  = loc[0];
    trace(i,1)  = loc[1];
    trace(i,2)  = sub_collapse_result["iv"];
    trace(i,3)  = sub_collapse_result["x_stat"];
    //monotony(c_stat) makes no sense in "All pairs" method
    if(mode == "J")
    {
      trace(i,4)  = sub_collapse_result["c_stat"];
      //trace(i,5)  =
    }
    trace(i,6)  = sub_collapse_result["ll"];
    trace(i,7)  = -2 * (trace(i,6) - ll0);
    trace(i,8)  = sub_collapse_result["lo_zscore"];
    sub_collapse_result = sub_collapse(freqMatrix_, method, mode);
  }

  return(trace);
}
