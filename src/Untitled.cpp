#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
int combine(NumericMatrix freqMatrix, int i, int j)
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
  freqMatrix = freqMatrix_;
  return(0);
}
