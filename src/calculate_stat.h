#ifndef _CALSTAT_h_
#define _CALSTAT_h_

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

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

inline double cal_iv(NumericVector good, NumericVector bad)
{
  double iv = 0;
  size_t n = good.size();
  for(size_t i = 0; i < n; ++i)
  {
    iv += (good[i] - bad[i]) * log(good[i] / bad[i]);
  }
  return(iv);
}

inline double cal_ll(NumericVector good, NumericVector bad)
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

inline double cal_log_odds_ratio_zscore(NumericVector good, NumericVector bad, int i, int j)
{
  double z_score;
  z_score = log((good[i] * bad[j]) / (good[j] * bad[i])) / sqrt(1 / good[i] + 1 / good[j] + 1 / bad[i] + 1 / bad[j]);
  return(z_score);
}

NumericMatrix combine(NumericMatrix freqMatrix, int i, int j)
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

void combineLabels(StringVector names, int i, int j)
{
  names[i] += ",";
  names[i] += names[j];
  for(int n = j + 1; n < names.size(); n ++)
  {
    names[n - 1] = names[n];
  }
  names[names.size() - 1] = " ";
}

double delta(NumericVector good, NumericVector bad, int i, int j, String method = "iv")
{
  double delta_ = 9999.0;
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
    double lo_zscore = cal_log_odds_ratio_zscore(good, bad, i, j);
    //under the condition that the two levels can be collapsed
    //we find the min adjust_lift(x_stat vs c_stat)
    if(abs(lo_zscore) > 1.96) //p_value <= 0.05
      return(9999.0);
    int n = good.size() - 1;
    NumericVector new_good(n);
    NumericVector new_bad(n);
    for(int k = 0; k < n; ++k)
    {
      if(k == i)
      {
        new_good[k] = good[i] + good[j];
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
    double c_stat = cal_c_stat(new_good, new_bad);
    delta_ = (cal_x_stat(new_good, new_bad) - c_stat) / c_stat;
    //if x_stat == c_stat then use vi maximum method
    if(delta_ <= 10e-6)
      return(9999.0);
  }

  return(delta_);
}

double delta_zerocell(NumericVector good, NumericVector bad, int zerolevel, int nonzerolevel)
{
  double sum_good = sum(good);
  double sum_bad  = sum(bad);
  double delta_   = (good[nonzerolevel] / sum_good - bad[nonzerolevel] / sum_bad) * (log(good[nonzerolevel] / sum_good) - log(bad[nonzerolevel] / sum_bad)) -
                    ((good[nonzerolevel] + good[zerolevel]) / sum_good - (bad[nonzerolevel] + bad[zerolevel]) / sum_bad) * (log((good[nonzerolevel] + good[zerolevel]) / sum_good) - log((bad[nonzerolevel] + bad[zerolevel]) / sum_bad));
  return(delta_);
}

#endif
