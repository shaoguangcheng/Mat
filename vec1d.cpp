
template <class T>
vec1d<T>::vec1d() 
  : nElem(0), refCount(new int(1)), data(NULL)
{}

template <class T>
vec1d<T>::vecld(int n, const T* data_)
  : nElem(n), refCount(new int(1))
{
  if(n <= 0){
    DEBUGMSG("invalid input");
    exit(-1);
  }

  data = new T [nElem];
  assert(data != NULL);
  
  if(NULL == data_)
    memset(data, 0, sizeof(T)*nElem);
  else
    memcpy(data, data_, nElem);
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
