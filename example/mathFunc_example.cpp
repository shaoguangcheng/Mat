#include "../src/mathFunc.hpp"

void testVec()
{
  std::cout << "=================" << std::endl;

  vec1d<float> v1(5);
  
  std::cout << "v1: " << v1;
  std::cout << "ref: " << v1.ref() << std::endl;

  vec1d<float> v2(5, 3);

  std::cout << "v2: " << v2;
  std::cout << "ref: " << v2.ref() << std::endl;

  v2 = v1;
  std::cout << "v1: " << v1;
  std::cout << "v2: " << v2;
  std::cout << "v1 ref: " << v1.ref() << std::endl;
  std::cout << "v2 ref: " << v2.ref() << std::endl;

  vec1d<float> v3(v2);
  std::cout << "v2: " << v2;
  std::cout << "v3: " << v3;
  std::cout << "v2 ref: " << v2.ref() << std::endl;
  std::cout << "v3 ref: " << v3.ref() << std::endl;  
}

void testScaleVec()
{
  std::cout << "=====================" << std:: endl;

  vec1d<float> v1(5);

  for(int i = 0; i < v1.size(); ++i)
    v1(i) = float(i+1);
  
  std::cout << "Before scale : " << v1 << std::endl;
  
  Mat_Scale(v1, float(2.0));

  std::cout << "After scale  : " << v1 << std::endl;

  vec1d<float> v2(v1.size(), 2), v3(v1.size());
  
  Mat_Mul(v1, v2, v3);

  float max;
  int index;
  Mat_Max(v3, max, index);
  std::cout << max << std::endl;
}


int main()
{
  testVec();
  testScaleVec();

  return 0;
}
