#include <Rcpp.h>
#include "calculate_stat.h"
using    namespace Rcpp;
using    namespace std;

double binary_split(NumericMatrix freqMatrix);

//[[Rcpp::export]]
NumericMatrix Collapse(NumericMatrix freqMatrix, NumericMatrix trace, int row_indx = 1, String method = "max_iv", String mode = "J")
{
  NumericVector good = freqMatrix.column(0);
  NumericVector bad  = freqMatrix.column(1);
  size_t n           = good.size();
  double min_delta   = R_PosInf, delta_;
  int left, right;

  //calculate initial log-likehood and max-binary-informance-value
  if(row_indx == 1)
  {
    trace(0, 2)  = min(good + bad);
    trace(0, 8)  = cal_ll(good, bad);
    trace(0, 9)  = 1;
    trace(0, 11) = mode == "J" ? binary_split(freqMatrix):NA_REAL;
    NumericVector good__ = good / sum(good);
    NumericVector bad__  = bad / sum(bad);
    trace(0, 3)  = cal_iv(good__, bad__);
    trace(0, 4)  = 0;
    trace(0, 6)  = mode == "J" ? cal_c_stat(good__, bad__):NA_REAL;
    trace(0, 5)  = n == 2 ? trace(0, 6):cal_x_stat(good__, bad__);
    trace(0, 7)  = mode == "J" ? (n == 2 ? 0:((trace(0, 5) - trace(0, 6)) / (trace(0, 6) * (n - 2)))):NA_REAL;
  }

  //find then best collapse
  if(mode == "J")
  {
    for(size_t i = 0; i < n - 1; ++i)
    {
      delta_ = delta(good, bad, i, i + 1, method);
      //when using method 'linear',if it's already linear
      //or there's no levels can be collapsed unsignificantly, use method 'max_iv'
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

  if(min_delta == 9999.0) return(Collapse(freqMatrix, trace, row_indx, "max_iv", "J"));

  //get new freqmatrix
  NumericMatrix freqMatrix_ = combine(freqMatrix, left, right);
  NumericVector new_good    = freqMatrix_.column(0);
  NumericVector new_bad     = freqMatrix_.column(1);

  //trace
  trace(row_indx, 0)  = left;
  trace(row_indx, 1)  = right;
  trace(row_indx, 2)  = min(new_good + new_bad);
  trace(row_indx, 10) = cal_log_odds_ratio_zscore(good, bad, left, right);
  trace(row_indx, 8)  = cal_ll(new_good, new_bad);
  new_good            = new_good / sum(new_good);
  new_bad             = new_bad / sum(new_bad);
  trace(row_indx, 3)  = n ==2 ? 0:cal_iv(new_good, new_bad);
  trace(row_indx, 4)  = (trace(row_indx, 3) - trace(0, 3)) * 100 / trace(0, 3);
  trace(row_indx, 6)  = mode == "J" ? cal_c_stat(new_good, new_bad):NA_REAL;
  trace(row_indx, 5)  = n == 2 ? trace(row_indx, 6):cal_x_stat(new_good, new_bad);
  trace(row_indx, 7)  = mode == "J" ? (n == 2 ? 0:((trace(row_indx, 5) - trace(row_indx, 6)) / (trace(row_indx, 6) * (n - 2)))):NA_REAL;
  trace(row_indx, 12) = method == "max_iv" ? 1:(method == "max_likehood" ? 2:3);

  if(n == 2) return(trace);
  return(Collapse(freqMatrix_, trace, ++ row_indx, method, mode));
}


double binary_split(NumericMatrix freqMatrix)
{
  NumericVector good = wrap(freqMatrix.column(0));
  NumericVector bad  = wrap(freqMatrix.column(1));
  size_t nr          = freqMatrix.nrow();
  double max_iv      = R_NegInf, iv;

  for(size_t i = 1; i <= nr - 1; ++i)
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

//[[Rcpp::export]]
NumericMatrix CollapseZeroCells(NumericMatrix freqMatrix, IntegerMatrix trace, String mode, int iter = 0)
{
  int nrow = freqMatrix.rows();

  int zeroflag = 0;
  int row, col;
  for(row = 0; row < nrow; row ++)
  {
    for(col = 0; col < 2; col ++)
    {
      if(freqMatrix(row, col) == 0)
      {
        zeroflag ++;
        break;
      }
    }
    if(zeroflag > 0) break;
  }

  if(zeroflag == 0 or nrow < 2 or iter > 100) return(freqMatrix);

  NumericVector good = freqMatrix.column(0);
  NumericVector bad  = freqMatrix.column(1);

  int left = -1, right = -1;     //calculate the left and right which will be collapsed
  if(mode == "J")
  {
    if(row == 0) left = row;    //the first row
    else if(row == nrow - 1) left = row - 1;    //the last row
    else
    {
      if(freqMatrix(row + 1, col) == 0)    //if the next row also has zero in the same column
      {                                    //collapse them directly, because can not use "max_iv","max_likehood" methods
        left = row;
      }
      else if(freqMatrix(row + 1, 1 - col) == 0)    //if the next row has zero in the different column
      {                                             //collapse the front row directly
        left = row - 1;
      }
      else
      {
        double delta_1 = delta_zerocell(good, bad, row, row - 1);
        double delta_2 = delta_zerocell(good, bad, row, row + 1);
        left = delta_1 < delta_2 ? (row - 1):row;
      }
    }
    right = left + 1;
  }
  else
  {
    double min_delta = R_PosInf;
    left = row;
    for(int row_ = 0; row_ < nrow; row_ ++)
    {
      if(row_ != row)
      {
        if(freqMatrix(row_, col) == 0)
        {
          right = row_;    //preferentialy collapse the levels which have zero in the same column
          break;
        }
        else if(freqMatrix(row_, 1 - col) > 0)  //make sure not to collapse the levels which has zero in different column
        {
          double delta_ = delta_zerocell(good, bad, row, row_);
          if(delta_ < min_delta)
          {
            min_delta = delta_;
            right = row_;
          }
        }
      }
    }
    if(right == -1) return(freqMatrix);      //in mode "A", there will be two levels and the diagnol elements are zeros
  }


  NumericMatrix freqMatrix_ = combine(freqMatrix, left, right);
  List dimnames = freqMatrix.attr("dimnames");
  StringVector labels = dimnames[0];
  combineLabels(labels, left, right);
  StringVector labels_(freqMatrix_.nrow());
  for(int i = 0; i < freqMatrix_.nrow(); ++ i) labels_[i] = labels[i];
  freqMatrix_.attr("dimnames") = List::create(labels_, dimnames[1]);
  iter ++;
  return(CollapseZeroCells(freqMatrix_, trace, mode, iter));
}


//[[Rcpp::export]]
NumericMatrix combineResults(NumericMatrix freqMatrix, NumericVector left, NumericVector right)
{
  NumericMatrix freqMatrix_ = clone(freqMatrix);
  List dimnames             = freqMatrix_.attr("dimnames");
  StringVector labels       = dimnames[0];
  if(labels.size() == 1) return freqMatrix_;
  for(int i = 0; i < left.size(); ++i)
  {
    combineLabels(labels, left[i], right[i]);
    freqMatrix_ = combine(freqMatrix_, left[i], right[i]);
  }

  StringVector labels_(freqMatrix.nrow() - left.size());
  for(int i = 0; i < labels_.size(); i ++) labels_[i] = labels[i];
  freqMatrix_.attr("dimnames") = List::create(labels_, dimnames[1]);
  return freqMatrix_;
}

//[[Rcpp::export]]
double linearity(NumericVector good, NumericVector bad)
{
  double c_stat    = cal_c_stat(good, bad);
  double x_stat    = cal_x_stat(good, bad);
  double adj_lift  = (x_stat - c_stat) / (c_stat * (good.size() - 2));
  return adj_lift;
}
