#include "../src/mathFunc.hpp"

void testScaleVec()
{
  vec1d<float> v1(5);

  for(float i = 0; i < v1.size(); ++i)
    v1(i) = float(i+1);
  
  std::cout << "Before scale : " << v1 << std::endl;
  
  scale(v1, float(2.0));

  std::cout << "After scale  : " << v1 << std::endl;

  vec1d<float> v2(v1.size());
  Mat_Exp(v1, v2);
  v2.print();
}


int main()
{
  testScaleVec();

  return 0;
}
