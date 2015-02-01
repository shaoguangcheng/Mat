#include "../src/mathFunc.hpp"

void testScaleVec()
{
  vec1d<float> v1(5);

  for(int i = 0; i < v1.size(); ++i)
    v1(i) = float(i+1);
  
  std::cout << "Before scale : " << v1 << std::endl;
  
  Mat_Scale(v1, float(2.0));

  std::cout << "After scale  : " << v1 << std::endl;

  vec1d<float> v2(v1.size(), 2), v3(v1.size());

  Mat_Mul(v1, v2, v3);
  v3.print();
  
  std::cout << Mat_Sum(v3) << std::endl;
}


int main()
{
  testScaleVec();

  return 0;
}