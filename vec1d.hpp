#ifndef __VEC1D__
#define __VEC1D__

template <class T>
class vec1d;
template <class T>
ostream& operator << (ostream& out, const vec1d<T>& v) const;

template <class T>
class vec1d{
private:
  int nElem;
  int* refCount;

public:
  T* data;
  
public:
  vec1d();
  vec1d(int n, T* data_ = NULL);
  vec1d(const vec1d<T>& v);
  ~vec1d();
  
  vec1d<T>& operator = (const vec1d<T>& v);

  T& operator [] (int index);
  T& operator [] (int index) const;
  vecld<T> range(int start, int end) const;

  inline int size() const {return nElem;}
  inline int ref() const {return *refCount;}
  
  T sum() const;
  T mean() const;
  T L1Norm() const;
  T L2Norm() const;

  void addScalar(T s);
  void addVec1d(const vec1d<T>& v);
  T dot(const vec1d<T>& v) const; // dot(x, y) = x*transpos(y)
  
  //only for two 3-element vector 
  vec1d<T> cross(const vec1d<T>& v) const; 

  friend ostream& operator << (ostream& out, const vec1d<T>& v) const;
  void print() const;
};

#include "vec1d.cpp"

#endif // end of vec1d
