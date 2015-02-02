/*!
 * \file mat2d.hpp
 * \breif 2-dim matrix is defined and implemented in this file.
 *   The mat2d supports range index, slice operations.
 *   To manage the memory effectively, reference count trick is used here.
 *
 * \author Shaoguang Cheng. From Xi'an, China
 * \date   2/2/2015
 */

#ifndef __MAT2D__
#define __MAT2D__

#include "global.h"

template <class T>
class vec1d;
template <class T>
std::ostream& operator << (std::ostream& out, const mat2d<T>& m);

template <class T>
class mat2d{
private:
  int nCol;
  int nRow;

public:
  int* refCount;
  T* data;

public:
  mat2d();
  mat2d(int rows, int cols, T* data_ = NULL);
  mat2d(int rows, int cols, const T& val)
  mat2d(const mat2d<T>& m);
  ~mat2d();

  static mat2d<T> zeros(int rows, int cols);
  static mat2d<T> eye(int rows, int cols);

  mat2d<T>& operator = (const mat2d<T>& m);
  mat2d<T> deepCopy() const;

  inline T& operator () (int indexX, int indexY);
  inline const T& operator () (int indexX, int indexY) const;
  
  inline int rows() const;
  inline int cols() const;
  inline int size() const;
  inline int ref() const;

  void reshape(int rows, int cols);

  vec1d<T> row(int index) const;
  vec1d<T> col(int index) const;
  mat2d<T> rowRange(int start, int end, int inc = 1) const;
  mat2d<T> colRange(int start, int end, int inc = 1) const;
  mat2d<T> subMat(int leftUpperX, int leftUpperY, 
		  int rightUpperX, int rightUpperY) const;
  
  bool operator == (const mat2d<T>& m) const;
  
  friend std::ostream& operator << <T> (std::ostream& out, const mat2d<T>& v);
  void print(char sep = ' ') const;

private:
  void decreaseUse();
};


#include "mat2d.cpp"
#endif //end of __mat2d
