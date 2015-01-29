#ifndef __MATRIX2D__
#define __MATRIX2D__

template <class T>
class matrix2d{
private:
  int* ref;
  int nCol;
  int nRow;

public:
  T* data;

public:
  matrix2d();
  matrix2d(int cols, int rows, T* data_ = NULL);
  matrix2d(const matrix2d<T>& m);
  ~matrix2d();

  matrix2d<T>& operator = (const matrix2d<T>& m);
  
  inline int rows() const {return nRow};
  inline int cols() const {return nCol;}
  inline int total() const {return nRow*nCol;};
  
  
};

#endif //end of __MATRIX2D
