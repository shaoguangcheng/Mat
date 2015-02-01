
template <class T>
vec1d<T>::vec1d() 
  : nElem(0), refCount(new int(1)), data(NULL)
{}

template <class T>
vec1d<T>::vec1d(int n, const T* data_)
  : nElem(n), refCount(new int(1)), data(new T [nElem])
{
  if(n <= 0){
    DEBUGMSG("size must be positive");
    exit(-1);
  }

  assert(data != NULL);
  
  if(NULL == data_)
    memset(data, 0, sizeof(T)*nElem);
  else
    memcpy(data, data_, nElem*sizeof(T));
}

template <class T>
vec1d<T>::vec1d(int n, T val) 
  : nElem(n), refCount(new int(1)), data(new T[nElem])
{
  if(n <= 0){
    DEBUGMSG("size must be positive");
    exit(-1);
  }

  assert(data != NULL);
  for(int i = 0; i < nElem; ++i)
    data[i] = val;
}

template <class T>
vec1d<T>::~vec1d()
{
  decreaseUse();
}

template <class T>
vec1d<T>::vec1d(const vec1d<T>& v)
  : nElem(v.size()), refCount(v.refCount), data(v.data)
{
  ++(*refCount);
}

template <class T>
vec1d<T>& vec1d<T>::operator = (const vec1d<T>& v)
{
  if(this != &v){
    decreaseUse();
    ++(*(v.refCount));
    nElem = v.size();
    refCount = v.refCount;
    data = v.data;
  }

  return *this;
}

template <class T>
void vec1d<T>::decreaseUse()
{
  if(*refCount == 1){
    delete [] data;
    delete refCount;
    data = NULL;
    refCount = NULL;
  }
  else
    --(*refCount);
}

template <class T>
vec1d<T> vec1d<T>::deepCopy() const
{
  vec1d<T> v(nElem, data);
  return v;
}

template <class T>
T& vec1d<T>::operator [] (int index)
{
  if(index >= nElem || index < 0){
    DEBUGMSG("index out of range");
    exit(-1);
  }

  return data[index];
}

template <class T>
const T& vec1d<T>::operator [] (int index) const
{
  if(index >= nElem || index < 0){
    DEBUGMSG("index out of range");
    exit(-1);
  }
  
  return data[index];
}

template <class T>
T& vec1d<T>::operator () (int index)
{
  if(index >= nElem || index < 0){
    DEBUGMSG("index out of range");
    exit(-1);
  }

  return data[index];
}

template <class T>
const T& vec1d<T>::operator () (int index) const
{
  if(index >= nElem || index < 0){
    DEBUGMSG("index out of range");
    exit(-1);
  }
  
  return data[index];
}

template <class T>
vec1d<T> vec1d<T>::operator () (int start, int end) const
{
  if(end <= start){
    DEBUGMSG("end index must be larger than start index");
    exit(-1);
  }

  if(start < 0 || end >= nElem){
    DEBUGMSG("index out of range");
    exit(-1);
  }

  if(end == -1)
    end = nElem;
  
  vec1d<T> v(end - start + 1, data + start);
  
  return v;
}

template <class T>
vec1d<T> vec1d<T>::range(int start, int end) const
{
  return (*this)(start, end);
}

template <class T>
void vec1d<T>::print(char sep) const
{
  if(NULL == data){
    std::cout << "[ NULL ]" << std::endl; 
  }
  else{
    std::cout << "[ ";
    for(int i = 0; i < nElem; ++i)
      std::cout << data[i] << sep;
    std::cout << "]" <<std::endl;
  }
}

template <class T>
std::ostream& operator << (std::ostream& out, const vec1d<T>& v) 
{
  if(NULL == v.data){
    out << "[ NULL ]" << std::endl; 
  }
  else{
    out << "[ ";
    for(int i = 0; i < v.nElem; ++i)
      out << v.data[i] << " ";
    out << "]"<< std::endl;
  }

  return out;
}
