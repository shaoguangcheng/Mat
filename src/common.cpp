
template <class T>
bool equal(const T& x, const T& y)
{
  T diff = x - y;

  return diff > -1e-5 && diff < 1e-5; 
}
