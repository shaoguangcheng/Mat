//#include "mathFunc.hpp"

template <class T>
void scale(vec1d<T>& X, const T& alpha)
{
  int size = X.size();
  for(int i = 0; i < size; ++i)
    X[i] *= alpha;
}

template <>
void scale<float>(vec1d<float>& X, const float& alpha)
{
  cblas_sscal(X.size(), alpha, X.data, 1);
}

template <>
void scale<double> (vec1d<double>& X, const double& alpha)
{
  cblas_dscal(X.size(), alpha, X.data, 1);
}

