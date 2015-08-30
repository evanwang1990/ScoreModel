#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
double cal_c_stat(NumericVector good, NumericVector bad)
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
double cal_x_stat(NumericVector good, NumericVector bad)
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
NumericMatrix combine(NumericMatrix freqMatrix, int i, int j)
{
  i--; j--;
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
  double delta;
  if(method == "iv")
  {
    int sum_good = sum(good);
    int sum_bad = sum(bad);
    delta = (good[i] / sum_good - bad[i] / sum_bad) * (log(good[i] / sum_good) - log(bad[i] / sum_bad)) +
            (good[j] / sum_good - bad[j] / sum_bad) * (log(good[j] / sum_good) - log(bad[j] / sum_bad)) -
            ((good[i] + good[j]) / sum_good - (bad[i] + bad[j]) / sum_bad) * (log((good[i] + good[j]) / sum_good) - log((bad[i] + bad[j]) / sum_bad));
  }
  else
  {
    int new_good = good[i] + good[j];
    int new_bad = bad[i] + bad[j];
    delta = good[i] * log(good[i] / (good[i] + bad[i])) + bad[i] * log(bad[i] / (good[i] + bad[i])) +
            good[j] * log(good[j] / (good[j] + bad[j])) + bad[j] * log(bad[j] / (good[j] + bad[j])) -
            new_good * log(new_good / (new_good + new_bad)) + new_bad * log(new_bad / (new_good + new_bad));
  }

  return(delta);
}

//[[Rcpp::export]]
NumericVector collapse
