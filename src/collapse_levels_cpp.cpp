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

//[[Rcpp::export]]
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

//[[Rcpp::export]]
inline double cal_log_odds_ratio_zscore(NumericVector good, NumericVector bad, int i, int j)
{
  double z_score;
  z_score = log((good[i] * bad[j]) / (good[j] * bad[i])) / sqrt(1 / good[i] + 1 / good[j] + 1 / bad[i] + 1 / bad[j]);
  return(z_score);
}

//[[Rcpp::export]]
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
  names[i] += "+";
  names[i] += names[j];
  for(int n = j + 1; n < names.size(); n ++)
  {
    names[n - 1] = names[n];
  }
  names[names.size() - 1] = " ";
}

//[[Rcpp::export]]
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
    if(abs(lo_zscore) > 1.64) return(9999.0);
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
    if(delta_ <= 10e-6) return(9999.0);
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
      if(delta_ == 9999.0) return(sub_collapse(freqMatrix));
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
  double c_stat = method == "J" ? cal_c_stat(new_good, new_bad):NA_REAL;
  double adj_lift = method == "J" ? (x_stat - c_stat) / (c_stat * (n - 2)):NA_REAL;


  return(List::create(Named("freqMatrix") = freqMatrix_,
                      Named("stat") = NumericVector::create(loc[0], loc[1], iv, x_stat, c_stat, adj_lift, ll, NA_REAL, lo_zscore, NA_REAL, NA_REAL)));
}

//[[Rcpp::export]]
double binary_split(NumericMatrix freqMatrix)
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad = wrap(freqMatrix.column(1));
  size_t nr = freqMatrix.nrow();
  double max_iv = R_NegInf;
  double iv;
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
