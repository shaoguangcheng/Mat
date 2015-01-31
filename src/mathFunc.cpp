//#include "mathFunc.hpp"

//----------------------level 0----------------------

template <>
void Mat_Mul<float>(const vec1d<float>& X,
		    const float& alpha,
		    vec1d<float>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_scopy(size, X.data, 1, Y.data, 1);
  cblas_sscal(size, alpha, Y.data, 1);
}

template <>
void Mat_Mul<double>(const vec1d<double>& X, 
		     const double& alpha,
		     vec1d<double>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_dcopy(size, X.data, 1, Y.data, 1);
  cblas_dscal(size, alpha, Y.data, 1);
}

template <>
void Mat_Div<float>(const vec1d<float>& X,
		    const float& alpha,
		    vec1d<float>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_scopy(size, X.data, 1, Y.data, 1);
  cblas_sscal(size, float(1/alpha), Y.data, 1);
}

template <>
void Mat_Div<double>(const vec1d<double>& X, 
		     const double& alpha,
		     vec1d<double>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_dcopy(size, X.data, 1, Y.data, 1);
  cblas_dscal(size, double(1/alpha), Y.data, 1);
}


//----------------------level 1-----------------------

template <>
void Mat_Scale<float>(vec1d<float>& X, const float& alpha)
{
  cblas_sscal(X.size(), alpha, X.data, 1);
}

template <>
void Mat_Scale<double> (vec1d<double>& X, const double& alpha)
{
  cblas_dscal(X.size(), alpha, X.data, 1);
}

template <class T>
void Mat_axpy(const vec1d<T>& X, const T& alpha, vec1d<T>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }
  
  vec1d<T> Z(size);

  Mat_Mul(X, alpha, Z);
  Mat_Add(Z, Y, Y);
}

template <>
void Mat_axpy<float>(const vec1d<float>& X, const float& alpha, vec1d<float>& Y)
{
  cblas_saxpy(X.size(), alpha, X.data, 1, Y.data, 1);
}

template <>
void Mat_axpy<double>(const vec1d<double>& X, const double& alpha, vec1d<double>& Y)
{
  cblas_daxpy(X.size(), alpha, X.data, 1, Y.data, 1);
}
