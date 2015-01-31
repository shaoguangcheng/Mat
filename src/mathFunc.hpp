#ifndef __MATH_FUNC__
#define __MATH_FUNC__

#include "vec1d.hpp"

#include <cblas.h>
#include <math.h>

// used to represent transpose operations on a matrix
enum MAT_TRANSPOSE{
  MatNoTrans = 111, // for compatible with cblas
  MatTrans = 112,
  MatConjTrans = 113
};

// used to indicate the order of mat-mat multiply
enum MAT_SIDE{
  MatLeft = 141,
  MatRight = 142
};

// some element-wise operations
// such as Y[i] = exp(X[i])
#define MAT_UNARY_FUNC(name, operation)	 \
  template <class Dtype>			 \
  void Mat_##name(const Dtype& X,		 \
		  Dtype& Y)			 \
  {						 \
    if(X.size() != Y.size()){			 \
      DEBUGMSG("length must be same");		 \
      exit(-1);					 \
    }						 \
						 \
    int size = X.size(); 			 \
    for(int i = 0; i < size; ++i){		 \
      operation;				 \
    }						 \
  }						 

MAT_UNARY_FUNC(Exp, Y[i] = exp(X[i]));
MAT_UNARY_FUNC(Sqrt, Y[i] = sqrt(X[i]));
MAT_UNARY_FUNC(Abs, Y[i] = fabs(X[i]));
MAT_UNARY_FUNC(Cos, Y[i] = cos(X[i]));
MAT_UNARY_FUNC(Sin, Y[i] = sin(X[i]));

#define MAT_BINARY_FUNC(name, operation)		\
  template <class Dtype>				\
  void Mat_##name(const Dtype& X1,			\
		  const Dtype& X2,			\
		  Dtype& Y)				\
  {							\
  int size = X1.size();					\
  assert(size == X2.size()				\
	 && size == Y.size());				\
  for(int i = 0; i < size; ++i){			\
    operation;						\
  }							\
}

MAT_BINARY_FUNC(Add, Y[i] = X1[i] + X2[i]);
MAT_BINARY_FUNC(Sub, Y[i] = X1[i] - X2[i]);
MAT_BINARY_FUNC(Mul, Y[i] = X1[i] * X2[i]);
MAT_BINARY_FUNC(Div, Y[i] = X1[i] / X2[i]);

//-----------------level 1------------------------

// X = X + alpha
template <class T>
void addScalar(vec1d<T>&X, const T& alpha);

// X = alpha*X
template <class T>
void scale(vec1d<T>& X, const T& alpha);

// Y = alpha*X + Y
template <class T>
void Mat_axpy(const vec1d<T>& X, const T& alpha, vec1d<T>& Y);

// Y = X^2
template <class T>
vec1d<T> Mat_square(const vec1d<T>& X);

// sum(X)
template <class T>
T Mat_sum(const vec1d<T>& X);

// mean(X)
template <class T>
T Mat_mean(const vec1d<T>& X);



#include "mathFunc.cpp"
#endif // end of math_func
