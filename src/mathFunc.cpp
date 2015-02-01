
//----------------------level 0----------------------

// Y = X*X
template <class Dtype>
void Mat_Square(const Dtype& X, 
		Dtype& Y)
{
  int size = X.size();
  if(size != Y.size()){
    DEBUGMSG("length must be the same");
    exit(-1);
  } 

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x,
			  decltype(Y[0]) y){
      y = x * x;
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

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x,
			  decltype(Y[0]) y){
      y =  exp(x);
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

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x, 
			  decltype(Y[0]) y){
      y = sqrt(x);
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

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x,
			  decltype(Y[0]) y){
      y = fabs(x);
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

  MAT_UNARY_FUNC(X, Y, [](decltype(X[0]) x,
			  decltype(Y[0]) y){
      y = cos(x);
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
				decltype(X2[0]) x2,
				decltype(Y[0]) y){
		    y = x1+x2;
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
				decltype(X2[0]) x2,
				decltype(Y[0]) y){
		    y = x1-x2;
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
				decltype(X2[0]) x2,
				decltype(Y[0]) y){
		    y = x1/x2;
		  });
}

// Y = X + alpha (vector)
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

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](const Dtype& x, 
					const Dtype& a,
					Dtype& y){
      y = x + a;
    });
}

// Y = X - alpha (vector)
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

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](const Dtype& x, 
					const Dtype& a,
					Dtype& y){
      y = x - a;
    });
}

// Y = X * alpha (vector)
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

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](const Dtype& x, 
					const Dtype& a,
					Dtype& y){
      y = x * a;
    });
}


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

// Y = X / alpha (vector)
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

  MAT_BINARY_FUNC_PARAM(X, alpha, Y, [](const Dtype& x, 
					const Dtype& a,
					Dtype& y){
      y = x / a;
    });
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

// Y = X + alpha (matrix)
template <class Dtype>
void Mat_Add(const mat2d<Dtype>& X,
	     const Dtype& alpha,
	     mat2d<Dtype>& Y)
{
}

// Y = X - alpha (matrix)
template <class Dtype>
void Mat_Sub(const mat2d<Dtype>& X,
	     const Dtype& alpha,
	     mat2d<Dtype>& Y)
{
}

// Y = X * alpha (matrix)
template <class Dtype>
void Mat_Mul(const mat2d<Dtype>& X,
	     const Dtype& alpha,
	     mat2d<Dtype>& Y)
{
}

// Y = X / alpha (matrix)
template <class Dtype>
void Mat_Div(const mat2d<Dtype>& X,
	     const Dtype& alpha,
	     mat2d<Dtype>& Y)
{
}

// X = X + alpha (vector)
template <class Dtype>
void Mat_Inc(vec1d<Dtype>& X,
	     const Dtype& alpha)
{
  MAT_UNARY_FUNC(X, alpha, [](Dtype& x,
			      const Dtype& alpha){
		   x += alpha;
		 });
}


// X = X * alpha
template <class Dtype>
void Mat_Scale(vec1d<Dtype>& X,
	     const Dtype& alpha)
{
  MAT_UNARY_FUNC(X, alpha, [](Dtype& x,
			      const Dtype& alpha){
		   x *= alpha;
		 });
}

// X = X + alpha
template <class Dtype>
void Mat_Inc(mat2d<Dtype>& X,
	     const Dtype& alpha)
{
  MAT_UNARY_FUNC(X, alpha, [](Dtype& x,
			      const Dtype& alpha){
		   x += alpha;
		 });
}

// X = X * alpha
template <class Dtype>
void Mat_Scale(mat2d<Dtype>& X,
	     const Dtype& alpha)
{
  MAT_UNARY_FUNC(X, alpha, [](Dtype& x,
			      const Dtype& alpha){
		   x *= alpha;
		 });
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

#pragma omp parallel for reduction(+:s) if(size > 300)  
  for(int i = 0; i < size; ++i)
    s += X[i];

  return s;
}

template <class T>
T Mat_Mean(const vec1d<T>& X)
{
  return Mat_Sum(X)/X.size();
}

