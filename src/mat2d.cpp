template < class T>
mat2d<T>::mat2d()
  : nCol(0), nRow(0), nElem(nCol*nRow), data(NULL)
{} 

template <class T>
mat2d<T>::mat2d(int rows, int cols, T*data_)
  : nCol(cols), nRow(rows), nElem(nCol*nRow), data(new T [nElem])
{
  if(nCol <= 0 || nRow <= 0){
    DEBUGMSG("cols and rows must be positive");
    exit(-1);
  }

  assert(data != NULL);
  
  if(NULL == data_)
    memset(data, 0, sizeof(T)*nElem);
  else
    memcpy(data, data_, sizeof(T)*nElem);
}

template <class T>
mat2d<T>::mat2d(int rows, int cols, const T& val)
  : nCol(cols), nRow(rows), nElem(nCol*nRow), data(new T [nElem])
{
  if(rows <= 0 || cols <= 0){
    DEBUGMSG("cols and rows must be positive");\
    exit(-1);
  }

  assert(data != NULL);

  for(int i = 0; i < nElem; ++i)
    data[i] = val;
}

template <class T>
mat2d<T>::mat2d(const mat2d<T>& m)
  : nCol(m.nCol), nRow(m.nRow), nElem(m.nElem), use(m.use), data(m.data)
{
}

template <class T>
mat2d<T>::~mat2d()
{
  if(use.only())
    delete [] data;
}

template <class T>
mat2d<T>& mat2d<T>::operator = (const mat2d<T>& m)
{
  if(this != &m){
    if(use.reattach(m.use)){
      delete [] data;
      data = NULL;
    }
    
    nCol = m.cols();
    nRow = m.rows();
    nElem = m.size();
    data = m.data;
  }
  
  return *this;
}

template <class T>
mat2d<T> mat2d<T>::deepCopy() const
{
   return *(new mat2d<T>(nRow, nCol, data));
}

template <class T>
const T* & mat2d<T>::operator [] (int index) const
{
  return data + nCol*index;
}

template <class T>
T* & mat2d<T>::operator [] (int index)
{
  return data + nCol*index;
}

template <class T>
T& mat2d<T>::operator () (int index)
{
  return data[index]; // data(index/nCol, index%nCol)
}

template <class T>
const T& mat2d<T>::operator () (int index) const
{
  return data[index]; // data(index/nCol, index%nCol)
}

template <class T>
T& mat2d<T>::operator () (int indexX, int indexY)
{
  return data[indexX*nCol+indexY];
}

template <class T>
const T& mat2d<T>::operator () (int indexX, int indexY) const
{
  return data[indexX*nCol+indexY];
}

template <class T>
int mat2d<T>::rows() const
{
  return nRow;
}

template <class T>
int mat2d<T>::cols() const
{
  return nCol;
}

template <class T>
int mat2d<T>::size() const
{
  return nElem;
}

template <class T>
int mat2d<T>::ref() const
{
  return use.getCount();
}

template <class T>
void mat2d<T>::reshape(int rows, int cols)
{
  if(rows*cols != nElem){
    DEBUGMSG("size can not be changed while reshaping");
    exit(-1);
  }

  nCol = cols;
  nRow = rows;
}

template <class T>
vec1d<T> mat2d<T>::row(int index) const
{
  if(index < 0 || index >= nRow){
    DEBUGMSG("out of range");
    exit(-1);
  }

  vec1d<T> v(nCol, data+index*nCol);
  return v; 
}

template <class T>
vec1d<T> mat2d<T>::col(int index) const
{
  if(index < 0 || index >= nCol){
    DEBUGMSG("out of range");
    exit(-1);
  }

  vec1d<T> v(nRow);
  T* data_ = data + index;

  for(int i = 0; i < nRow; ++i){
    v[i] = *data_;
    data_ = data_ + nCol;
  }

  return v;
}

template <class T>
mat2d<T> mat2d<T>::rowRange(int start, int end) const
{
  int newRow = end - start + 1;

  if(start < 0 || end >= nRow ||
     newRow <= 0 || newRow > nRow){
    DEBUGMSG("Invalid input");
    exit(-1);
  }

  mat2d<T> m(newRow, nCol, data+nCol*start);

  return m;
}

template <class T>
mat2d<T> mat2d<T>::colRange(int start, int end) const
{
  int newCol = end - start + 1;
  
  if(start < 0 || end >= nCol || 
     newCol <= 0 || newCol > nCol){
    DEBUGMSG("Invalid input");
    exit(-1);
  }

  mat2d<T> m(nRow, newCol);
  T* newData_ = m.data;
  T* data_ = data + start;

  int tmp = 0, tmp_ = 0;
  for(int i = 0; i < nRow; ++i){
    tmp  = i*newCol;
    tmp_ = i*nCol; 
    for(int j = 0; j < newCol; ++j){
      newData_[tmp+j] = *(data_ + tmp_ + j);
    }
  }

  return m;
}

template <class T>
mat2d<T> mat2d<T>::subMat(int leftUpperX, int leftUpperY, 
			  int rightBottomX, int rightBottomY) const
{
  int newRow = rightBottomY - leftUpperY;
  int newCol = rightBottomX - leftUpperX;

  if(leftUpperX < 0 || leftUpperX > nRow ||
     leftUpperY < 0 || leftUpperY > nCol ||
     newRow <= 0 || newRow > nRow ||
     newCol <=0 || newCol > nCol){
    DEBUGMSG("Invalid input");
    exit(-1);
  }

  mat2d<T> m(newRow, newCol);
  T* newData = m.data;
  T* data_ = data + nCol*leftUpperX + leftUpperY;
  
  int tmp = 0, tmp_ = 0;
  for(int i = 0; i < newRow; ++i){
    tmp = i*newCol;
    tmp_ = i*nCol;
    for(int j = 0; j < newCol; ++j){
      newData[tmp+j] = *(data_ + tmp_ + j);
    }
  }

  return m;
}

template <class T>
bool mat2d<T>::operator == (const mat2d<T>& m) const
{
  if(m.rows() != nRow || m.cols() != nCol || m.size() != nElem)
    return false;
  
  for(int i = 0; i < nElem; ++i){
    if(data[i] != m(i))
      return false;
  }

  return true;
}

template <class T>
static mat2d<T> zeros(int rows, int cols)
{
  return *(new mat2d<T>(rows, cols, 0));
}

template <class T>
static mat2d<T> eye(int rows, int cols)
{
  mat2d<T> m(rows, cols);

  int tmp = rows > cols ? cols : rows;

  for(int i = 0; i < tmp; ++i)
    m(i, i) = static_cast<T>(1);

  return m;
}

template <class T>
void mat2d<T>::print(char sep) const
{
  if(nElem == 0)
    std::cout << "[NULL]" << std::endl;
  else{
    std::cout << "[ ";
    for(int i = 0; i < nElem; ++i)
      std::cout << data[i] << sep;
    std::cout << "]" << std::endl;
  }
}

template <class T>
std::ostream& operator << (std::ostream& out, const mat2d<T>& m) 
{
  if(NULL == m.data){
    out << "[ NULL ]" << std::endl; 
  }
  else{
    out << "[ ";
    for(int i = 0; i < m.nElem; ++i)
      out << m.data[i] << " ";
    out << "]"<< std::endl;
  }

  return out;
}

