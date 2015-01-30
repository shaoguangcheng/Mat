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

public:
  T* data;
  int* refCount;
  
public:
  vec1d();

  // create val with n elements and use an array to intialize 
  vec1d(int n, const T* data_ = NULL);

  // create val with n elements and initailize it with val
  vec1d(int n, T val);

  vec1d(const vec1d<T>& v);
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

  inline int size() const {return nElem;}
  inline int ref() const {return *refCount;}
  
  friend std::ostream& operator << <T>(std::ostream& out, const vec1d<T>& v);
  void print(char sep = ' ') const;

private:
  void decreaseUse();
};

#include "vec1d.cpp"

#endif // end of vec1d
