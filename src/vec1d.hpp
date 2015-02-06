/*!
 * \file vec1d.hpp
 * \breif 1-dim vector is defined and implemented in this file.
 *   The vec1d supports range index, slice operations.
 *   To manage the memory effectively, reference count trick is used here.
 *
 * \author Shaoguang Cheng. From Xi'an, China
 * \date   1/2/2015
 */


#ifndef __VEC1D__
#define __VEC1D__

#include "global.h"

template <class T>
class vec1d;
template <class T>
std::ostream& operator << (std::ostream& out, const vec1d<T>& v);

template <class T>
class vec1d{
private:
  int nElem;
  useCount use;

public:
  T* data;
  
public:
  vec1d();

  // create val with n elements and use an array to intialize 
  vec1d(int n, const T* data_ = NULL);

  // create val with n elements and initailize it with val
  vec1d(int n, const T& val);

  vec1d(const vec1d<T>& v);

  vec1d(const std::vector<T>& v);

  ~vec1d();
  
  vec1d<T>& operator = (const vec1d<T>& v);

  vec1d<T> deepCopy() const;

  // access vector element
  inline T& operator [] (int index);
  inline const T& operator [] (int index) const;
  inline T& operator () (int index);
  inline const T& operator () (int index) const;
  inline vec1d<T> operator () (int start, int end) const;
  inline vec1d<T> range(int start, int end) const;

  bool operator == (const vec1d<T>& v) const;
 
  std::vector<T> toVector()const;

  inline int size() const {return nElem;}
  inline int ref() const {return use.getCount();}
 
  friend std::ostream& operator << <T>(std::ostream& out, const vec1d<T>& v);
  void print(char sep = ' ') const;

};

#include "vec1d.cpp"

#endif // end of vec1d
