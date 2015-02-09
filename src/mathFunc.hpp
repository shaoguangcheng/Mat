/*!
 * \file mathFunc.hpp
 * \breif Some vector and matrix functions are implemented in this 
 * file. The low level operations are based on cblas due to the efficiency
 * reason. OpenMp is also used but not optimized yet.
 *
 * \author Shaoguang Cheng. From Xi'an, China
 * \date   1/2/2015
 */

#ifndef __MATH_FUNC__
#define __MATH_FUNC__

#include "vec1d.hpp"
#include "mat2d.hpp"

#include <cblas.h>
#include <math.h>

#define MAT(X) const enum CBLAS_##X
#define Mat(X) Cblas##X 

//----------------------level 0--------------------------

// Due to the  

#if 0
// some element-wise operations
// such as Y[i] = exp(X[i])
#define MAT_UNARY_FUNC(name, operation)		 \
  template <class Dtype>			 \
  void Mat_##name(const Dtype& X,		 \
		  Dtype& Y)			 \
  {						 \
    if(X.size() != Y.size()){			 \
      DEBUGMSG("length must be same");		 \
      exit(-1);					 \
    }						 \
						 \
    int size = X.size();			 \
    for(int i = 0; i < size; ++i){		 \
      operation;				 \
    }						 \
  }						 

MAT_UNARY_FUNC(Square, Y[i] = X[i] * X[i]); // Y = X^2
MAT_UNARY_FUNC(Exp, Y[i] = exp(X[i]));   // Y = e^X
MAT_UNARY_FUNC(Sqrt, Y[i] = sqrt(X[i])); // Y = sqrt(X)
MAT_UNARY_FUNC(Abs, Y[i] = fabs(X[i]));  // Y = |X|
MAT_UNARY_FUNC(Cos, Y[i] = cos(X[i]));   // Y = cos(X)
MAT_UNARY_FUNC(Sin, Y[i] = sin(X[i]));   // Y = sin(X)
MAT_UNARY_FUNC(Tan, Y[i] = tan(X[i]));   // Y = tan(X)
MAT_UNARY_FUNC(Log, Y[i] = log(X[i]));   // Y = log(X)
MAT_UNARY_FUNC(Log10, Y[i] = log10(X[i]));   // Y = log10(X)

#define MAT_BINARY_FUNC(name, operation)		\
  template <class Dtype>				\
  void Mat_##name(const Dtype& X1,			\
		  const Dtype& X2,			\
		  Dtype& Y)				\
  {							\
    int size = X1.size();				\
    if(size != X2.size() || size != Y.size()){		\
      DEBUGMSG("length must be the same");		\
      exit(-1);						\
    }							\
    							\
    for(int i = 0; i < size; ++i){			\
      operation;					\
    }							\
  }

MAT_BINARY_FUNC(Add, Y[i] = X1[i] + X2[i]); // Y = X1 + X2 
MAT_BINARY_FUNC(Sub, Y[i] = X1[i] - X2[i]); // Y = X1 - X2
MAT_BINARY_FUNC(Mul, Y[i] = X1[i] * X2[i]); // Y = X1 * X2
MAT_BINARY_FUNC(Div, Y[i] = X1[i] / X2[i]); // Y = X1 / X2


#define MAT_UNARY_FUNC_PARAM(name, operation) \
  template <class Dtype>			      \
  void Mat_##name(const vec1d<Dtype>& X,	      \
		  const Dtype& alpha,		      \
		  vec1d<Dtype>& Y)		      \
  {						      \
    int size = X.size();			      \
    if(size != Y.size()){			      \
      DEBUGMSG("length must be the same");	      \
      exit(-1);					      \
    }						      \
    						      \
    for(int i = 0; i < size; ++i){		      \
      operation;				      \
    }						      \
  }						      
						      
/*  template <class DType>			      \
  void Mat_##name(const mat2d<Dtype>& X,	      \
		  const Dtype& alpha,		      \
		  mat2d<Dtype>& Y)		      \
  {						      \
    int size = X.size();			      \
    if(size != Y.size()){			      \
      DEBUGMSG("length must be the same");	      \
      exit(-1);					      \
    }						      \
						      \
    for(int i = 0; i < size; ++i){		      \
      operation;				      \
    }						      \
    }						      
  */

MAT_UNARY_FUNC_PARAM(Add, Y[i] = X[i] + alpha); // Y = X + alpha
MAT_UNARY_FUNC_PARAM(Sub, Y[i] = X[i] - alpha); // Y = X - alpha
MAT_UNARY_FUNC_PARAM(Mul, Y[i] = X[i] * alpha); // Y = alpha * X
MAT_UNARY_FUNC_PARAM(Div, Y[i] = X[i] / alpha); // Y = X / alpha

#define MAT_UNARY_FUNC_SELF(name, operation) \
  template <class Dtype>				\
  void Mat_##name(vec1d<Dtype>& X, const Dtype& alpha)	\
  {							\
    int size = X.size();				\
    for(int i = 0; i < size; ++i){			\
      operation;					\
    }							\
  }							

/*
  template <class Dtype>				\
  void Mat_##name(mat2d<Dtype>& X, const Dtype& alpha)	\
  {							\
  int size = X.size();					\
  for(int i = 0; i < size; ++i){			\
    operation;						\
  }							\
}
*/

MAT_UNARY_FUNC_SELF(Inc, X[i] += alpha);   // X = X + alpha;
MAT_UNARY_FUNC_SELF(Scale, X[i] *= alpha); // X = alpha * X

#else // else of 0

// some element-wise operations 
// for vector and matrix
template <class Dtype, class Func>
void MAT_UNARY_FUNC(const Dtype& X, 
		    Dtype& Y,
		    Func op);

// for vector and matrix
template <class Dtype, class Func>
void MAT_BINARY_FUNC(const Dtype& X1, 
		     const Dtype& X2,
		     Dtype& Y,
		     Func op);

// for vector only
template <class Dtype, class Func>
void MAT_BINARY_FUNC_PARAM(const vec1d<Dtype>& X,
			   const Dtype& alpha,
			   vec1d<Dtype>& Y,
			   Func op);

// for matrix only
template <class Dtype, class Func>
void MAT_BINARY_FUNC_PARAM(const mat2d<Dtype>& X,
			   const Dtype& alpha,
			   mat2d<Dtype>& Y,
			   Func op);

// for vector only
template <class Dtype, class Func>
void MAT_UNARY_FUNC_SELF(vec1d<Dtype>& X,
			 const Dtype& alpha,
			 Func op);

// for matrix only
template <class Dtype, class Func>
void MAT_UNARY_FUNC_SELF(mat2d<Dtype>& X,
			 const Dtype& alpha,
			 Func op);

// Y = X*X (vector and matrix)
template <class Dtype>
void Mat_Square(const Dtype& X, Dtype& Y);

// Y = e^X (vector and matrix)
template <class Dtype>
void Mat_Exp(const Dtype& X, Dtype& Y);

// Y = sqrt(X) (vector and matrix)
template <class Dtype>
void Mat_Sqrt(const Dtype& X, Dtype& Y);

// Y = |X| (vector and matrix)
template <class Dtype>
void Mat_Abs(const Dtype& X, Dtype& Y);

// Y = cos(X) (vector and matrix)
template <class Dtype>
void Mat_Cos(const Dtype& X, Dtype& Y);

// Y = sin(X) (vector and matrix)
template <class Dtype>
void Mat_Sin(const Dtype& X, Dtype& Y);

// Y = tan(X) (vector and matrix)
template <class Dtype>
void Mat_Tan(const Dtype& X, Dtype& Y);

// Y = log(X) (vector and matrix)
template <class Dtype>
void Mat_Log(const Dtype& X, Dtype& Y);

// Y = log10(X) (vector and matrix)
template <class Dtype>
void Mat_Log10(const Dtype& X, Dtype& Y);

// Y = X1 + X2 (vector and matrix)
template <class Dtype>
void Mat_Add(const Dtype& X1, const Dtype& X2, Dtype& Y);

// Y = X1 - X2 (vector and matrix)
template <class Dtype>
void Mat_Sub(const Dtype& X1, const Dtype& X2, Dtype& Y);

// Y = X1 * X2 (vector and matrix)
template <class Dtype>
void Mat_Mul(const Dtype& X1, const Dtype& X2, Dtype& Y);

// Y = X1/X2 (vector and matrix)
template <class Dtype>
void Mat_Div(const Dtype& X1, const Dtype& X2, Dtype& Y);

// Y = X + alpha (vector)
template <class Dtype>
void Mat_Add(const vec1d<Dtype>& X, const Dtype& alpha, vec1d<Dtype>& Y);

// Y = X - alpha (vector)
template <class Dtype>
void Mat_Sub(const vec1d<Dtype>& X, const Dtype& alpha, vec1d<Dtype>& Y);

// Y = X * alpha (vector)
template <class Dtype>
void Mat_Mul(const vec1d<Dtype>& X, const Dtype& alpha, vec1d<Dtype>& Y);

// Y = X / alpha (vector)
template <class Dtype>
void Mat_Div(const vec1d<Dtype>& X, const Dtype& alpha, vec1d<Dtype>& Y);

// Y = X + alpha (matrix)
template <class Dtype>
void Mat_Add(const mat2d<Dtype>& X, const Dtype& alpha, mat2d<Dtype>& Y);

// Y = X - alpha (matrix)
template <class Dtype>
void Mat_Sub(const mat2d<Dtype>& X, const Dtype& alpha, mat2d<Dtype>& Y);

// Y = X * alpha (matrix)
template <class Dtype>
void Mat_Mul(const mat2d<Dtype>& X, const Dtype& alpha, mat2d<Dtype>& Y);

// Y = X / alpha (matrix)
template <class Dtype>
void Mat_Div(const mat2d<Dtype>& X, const Dtype& alpha, mat2d<Dtype>& Y);

// X = X + alpha (vector)
template <class Dtype>
void Mat_Inc(vec1d<Dtype>& X, const Dtype& alpha);

// X = X * alpha (vector)
template <class Dtype>
void Mat_Scale(vec1d<Dtype>& X, const Dtype& alpha);

// X = X + alpha (matrix)
template <class Dtype>
void Mat_Inc(mat2d<Dtype>& X, const Dtype& alpha);

// X = X * alpha (matrix)
template <class Dtype>
void Mat_Scale(mat2d<Dtype>& X, const Dtype& alpha);


#endif // end of 0

//-----------------level 1------------------------

// Y = alpha*X + Y
template <class T>
void Mat_axpy(const vec1d<T>& X, const T& alpha, vec1d<T>& Y);

// Y = alpha*X + beta*Y
template <class T>
void Mat_axpby(const vec1d<T>& X, const T& alpha, 
	       vec1d<T>& Y,const T& beta);

// y = X^t*Y
template <class T>
T Mat_Dot(const vec1d<T>& X, vec1d<T>& Y);

// Y = ||X||2
template <class T>
T Mat_Norm2(const vec1d<T>& X);

// sum(X)
template <class T>
T Mat_Sum(const vec1d<T>& X);

// mean(X)
template <class T>
T Mat_Mean(const vec1d<T>& X);

// max(X)
template <class T>
void Mat_Max(const vec1d<T>& X, T& m, int& index);

// min(X)
template <class T>
void Mat_Min(const vec1d<T>& X, T& m, int& index);

//---------------------level 2--------------------------

//---------------------level 3--------------------------

// Y = alpha*X + Y
template <class T>
void Mat_axpy(const mat2d<T>& X, const T& alpha, mat2d<T>& Y);

// Y = alpha*X + beta*Y
template <class T>
void Mat_axpby(const mat2d<T>& X, const T& alpha, 
	       mat2d<T>& Y,const T& beta);

// C = alpha*op(A)*op(B)+beta*C
// op(A) = {A, A.t, A.h}
template <class T>
void Mat_gemm(MAT(TRANSPOSE) transA,
	      MAT(TRANSPOSE) transB,
	      const mat2d<T>& A,
	      const mat2d<T>& B,
	      const T& alpha,
	      const T& beta,
	      mat2d<T>& C);

// C = alpha*A*B+beta*C
// or
// C = alpha*B*A+beta*C 
// where A = A.t
template <class T>
void Mat_symm(MAT(SIDE) sideA,
	      MAT(UPLO) uploA,
	      const mat2d<T>& A,
	      const mat2d<T>& B,
	      const T& alpha,
	      const T& beta,
	      mat2d<T>& C);

// C = alpha*A*A'+beta*C (noTrans) 
// or
// C = alpha*A'*A+beta*C (Trans)
// where C is n by n
// where c is an n-by-n symmetric matrix; 
// a is an n-by-k matrix, if trans = 'N'or'n', 
// a is a k-by-n matrix, if trans = 'T'or't','C'or'c'
template <class T>
void Mat_syrk(MAT(UPLO) uploA,
	      MAT(TRANSPOSE) transA,
	      const mat2d<T>& A,
	      const T& alpha,
	      const T& beta,
	      mat2d<T>& C);

// B = alpha*op(A)*B
// or
// B = alpha*B*op(A)
// where op(A) = {A, A', Ah}. B is m by n
// where b is an m-by-n general matrix, and a is triangular; 
// op(a) must be an m-by-m matrix, if side = 'L'or'l' 
// op(a) must be an n-by-n matrix, if side = 'R'or'r
template <class T>
void Mat_trmm(MAT(SIDE) sideA,
	      MAT(UPLO) uploA,
	      MAT(TRANSPOSE) transA,
	      const mat2d<T>& A,
	      const T& alpha,
	      mat2d<T>& B);

// B = alpha*op(A(-1))*B
// or
// B = alpha*B*op(A(-1))
// where op(A) = {A, A', Ah}
// op(a)*x = alpha*b  or  x*op(a) = alpha*b, 
// where x and b are m-by-n general matrices, and a is triangular; 
// op(a) must be an m-by-m matrix, if side = 'L'or'l' 
// op(a) must be an n-by-n matrix, if side = 'R'or'r'.
template <class T>
void Mat_trsm(MAT(SIDE) sideA,
	      MAT(UPLO) uploA,
	      MAT(TRANSPOSE) transA,
	      const mat2d<T>& A,
	      const T& alpha,
	      mat2d<T>& B);

#include "mathFunc.cpp"
#endif // end of math_func
