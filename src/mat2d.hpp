#ifndef __MAT2D__
#define __MAT2D__

template <class T>
class mat2d{
private:
  int* ref;
  int nCol;
  int nRow;

public:
  T* data;

public:
  mat2d();
  mat2d(int cols, int rows, T* data_ = NULL);
  mat2d(const mat2d<T>& m);
  ~mat2d();

  mat2d<T>& operator = (const mat2d<T>& m);
  
  inline int rows() const {return nRow;}
  inline int cols() const {return nCol;}
  inline int total() const {return nRow*nCol;};
  
  
};

#endif //end of __mat2d
