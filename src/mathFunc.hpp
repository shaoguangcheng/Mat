#ifndef __MATH_FUNC__
#define __MATH_FUNC__

#include "vec1d.hpp"
#include "mat2d.hpp"

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

//----------------------level 0--------------------------
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
		    Func op)
{
  int size = X.size();

#pragma omp parallel for shared(X, size) if(size > 300)
  for(int i = 0; i < size; ++i){
    Y[i] = op(X[i]);
  }
}

// for vector and matrix
template <class Dtype, class Func>
void MAT_BINARY_FUNC(const Dtype& X1, 
		     const Dtype& X2,
		     Dtype& Y,
		     Func op)
{
  int size = X1.size();

#pragma omp parallel for shared(X1, X2, size) if(size > 300)
  for(int i = 0; i < size; ++i){
    op(X1[i], X2[i], Y[i]);
  }
}

// for vector only
template <class Dtype, class Func>
void MAT_BINARY_FUNC_PARAM(const vec1d<Dtype>& X,
			   const Dtype& alpha,
			   vec1d<Dtype>& Y,
			   Func op)
{
  int size = X.size();
  
#pragma omp parallel for shared(X, alpha, size) if(size > 300)
  for(int i = 0; i < size; ++i){
    Y[i] = op(X[i], alpha);
  }
}

// Y = X*X
template <class Dtype>
void Mat_Square(const Dtype& X, Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return x*x;
    });
}

// Y = e^X
template <class Dtype>
void Mat_Exp(const Dtype& X, 
	     Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return exp(x);
    });
}

// Y = sqrt(X)
template <class Dtype>
void Mat_Sqrt(const Dtype& X, 
	      Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return sqrt(x);
    });
}

// Y = |X|
template <class Dtype>
void Mat_Abs(const Dtype& X, 
	     Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return fabs(x);
    });
}

// Y = cos(X)
template <class Dtype>
void Mat_Cos(const Dtype& X, 
	     Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return cos(x);
    });
}

// Y = sin(X)
template <class Dtype>
void Mat_Sin(const Dtype& X, 
	     Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x){
      return sin(x);
    });
}

// Y = X1 + X2
template <class Dtype>
void Mat_Add(const Dtype& X1, 
	     const Dtype& X2, 
	     Dtype& Y)
{
  int size = X1.size();					
  if(size != X2.size() || size != Y.size()){		
    DEBUGMSG("length must be the same");		
    exit(-1);						
  }

  MAT_BINARY_FUNC(X1, X2, Y, [](decltype(X1[0]) x1,
				decltype(X2[0]) x2){
		    return x1+x2;
		  });
}

// Y = X1 - X2
template <class Dtype>
void Mat_Sub(const Dtype& X1, 
	     const Dtype& X2, 
	     Dtype& Y)
{
  int size = X1.size();					
  if(size != X2.size() || size != Y.size()){		
    DEBUGMSG("length must be the same");		
    exit(-1);						
  }

  MAT_BINARY_FUNC(X1, X2, Y, [](decltype(X1[0]) x1,
				decltype(X2[0]) x2){
		    return x1-x2;
		  });
}

// Y = X1 * X2
template <class Dtype>
void Mat_Mul(const Dtype& X1, 
	     const Dtype& X2, 
	     Dtype& Y)
{
  int size = X1.size();					
  if(size != X2.size() || size != Y.size()){		
    DEBUGMSG("length must be the same");		
    exit(-1);						
  }

  MAT_BINARY_FUNC(X1, X2, Y, []( decltype(X1[0]) x1,
				 decltype(X2[0]) x2,
				decltype(*Y.data) y){
		    y = x1*x2;
		  });
}

// Y = X1/X2 
template <class Dtype>
void Mat_Div(const Dtype& X1, 
	     const Dtype& X2, 
	     Dtype& Y)
{
  int size = X1.size();					
  if(size != X2.size() || size != Y.size()){		
    DEBUGMSG("length must be the same");		
    exit(-1);						
  }

  MAT_BINARY_FUNC(X1, X2, Y, [](decltype(X1[0]) x1,
				decltype(X2[0]) x2){
		    return x1/x2;
		  });
}

// Y = X + alpha
template <class Dtype>
void Mat_Add(const vec1d<Dtype>& X,
	     const Dtype& alpha,
	     vec1d<Dtype>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](Dtype x, Dtype a){
      return x + a;
    });
}

// Y = X - alpha
template <class Dtype>
void Mat_Sub(const vec1d<Dtype>& X,
	     const Dtype& alpha,
	     vec1d<Dtype>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](Dtype x, Dtype a){
      return x - a;
    });
}

// Y = X * alpha
template <class Dtype>
void Mat_Mul(const vec1d<Dtype>& X,
	     const Dtype& alpha,
	     vec1d<Dtype>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](Dtype x, Dtype a){
      return x * a;
    });
}

// Y = X / alpha
template <class Dtype>
void Mat_Div(const vec1d<Dtype>& X,
	     const Dtype& alpha,
	     vec1d<Dtype>& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  }

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](Dtype x, Dtype a){
      return x / a;
    });
}



#endif // end of 0
//-----------------level 1------------------------

// Y = alpha*X + Y
template <class T>
void Mat_axpy(const vec1d<T>& X, 
	      const T& alpha, 
	      vec1d<T>& Y);

// Y = alpha*X + beta*Y
template <class T>
void Mat_axpby(const vec1d<T>& X, 
	       const T& alpha, 
	       vec1d<T>& Y,
	       const T& beta);

// y = X^t*Y
template <class T>
T Mat_Dot(const vec1d<T>& X, 
	  vec1d<T>& Y);

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
void Mat_Max(const vec1d<T>& X, 
	     T& m, 
	     int& index = -1);

// min(X)
template <class T>
void Mat_Min(const vec1d<T>& X, 
	     T& m, 
	     int& index = -1);

#include "mathFunc.cpp"
#endif // end of math_func
