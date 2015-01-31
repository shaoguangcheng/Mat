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
void Mat_Scale<float>(vec1d<float>& X, 
		      const float& alpha)
{
  cblas_sscal(X.size(), alpha, X.data, 1);
}

template <>
void Mat_Scale<double> (vec1d<double>& X, 
			const double& alpha)
{
  cblas_dscal(X.size(), alpha, X.data, 1);
}

template <class T>
void Mat_axpy(const vec1d<T>& X, 
	      const T& alpha, vec1d<T>& Y)
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
void Mat_axpy<float>(const vec1d<float>& X, 
		     const float& alpha, 
		     vec1d<float>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_saxpy(size, alpha, X.data, 1, Y.data, 1);
}

template <>
void Mat_axpy<double>(const vec1d<double>& X, const double& alpha, vec1d<double>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_daxpy(size, alpha, X.data, 1, Y.data, 1);
}

template <class T>
void Mat_axpby(const vec1d<T>& X, 
	       const T& alpha, 
	       vec1d<T>& Y, 
	       const T& beta)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }
  
  Mat_Scale(Y, beta);
  Mat_axpy(X, alpha, Y);
}

template <>
void Mat_axpby<float>(const vec1d<float>& X,
		      const float& alpha,
		      vec1d<float>& Y,
		      const float& beta)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_sscal(size, beta, Y.data, 1);
  cblas_saxpy(size, alpha, X.data, 1, Y.data, 1);
}

template <>
void Mat_axpby<double>(const vec1d<double>& X,
		      const double& alpha,
		      vec1d<double>& Y,
		      const double& beta)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  cblas_dscal(size, beta, Y.data, 1);
  cblas_daxpy(size, alpha, X.data, 1, Y.data, 1);
}

template <class T>
T Mat_Dot(const vec1d<T>& X, vec1d<T>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  Mat_Mul(X, Y, Y);
  return Mat_Sum(Y);
}

template <>
float Mat_Dot<float>(const vec1d<float>& X, vec1d<float>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  return cblas_sdot(size, X.data, 1, Y.data, 1);
}

template <>
double Mat_Dot<double>(const vec1d<double>& X, vec1d<double>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }
  
  return cblas_ddot(size, X.data, 1, Y.data, 1);
}

template <class T>
T Mat_Norm2(const vec1d<T>& X)
{
  T s = 0;
  int size = X.size();

  vec1d<T> Y(size);

  Mat_Mul(X, X, Y);
  return sqrt(Mat_Sum(Y));
}

template <>
float Mat_Norm2<float>(const vec1d<float>& X)
{
  return cblas_snrm2(X.size(), X.data, 1);
}

template <>
double Mat_Norm2<double>(const vec1d<double>& X)
{
  return cblas_dnrm2(X.size(), X.data, 1);
}


template <class T>
T Mat_Sum(const vec1d<T>& X)
{
  T s = T(0);
  int size = X.size();
  
  for(int i = 0; i < size; ++i)
    s += X[i];

  return s;
}
