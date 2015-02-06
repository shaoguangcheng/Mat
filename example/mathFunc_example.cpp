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
 }

void testMat2d()
{
  float array [48];

  for(int i = 0; i < 48; ++i)
    array[i] = i+1;

  mat2d<float> m1(6, 8, array);

  std::cout << "m1 : "  << m1;
  std::cout << "ref: " << m1.ref() << std::endl;

  mat2d<float> m2(6, 8, 3);
  
  std::cout << "m2 : "  << m2;
  std::cout << "ref: " << m2.ref() << std::endl;
  
  m2 = m1;
  std::cout << "m1 : "  << m1;
  std::cout << "m2 : "  << m2;
  std::cout << "ref: " << m1.ref() << std::endl;
  std::cout << "ref: " << m2.ref() << std::endl;

  mat2d<float> m3(m1);

  std::cout << "m1 : "  << m1;
  std::cout << "m3 : "  << m3;
  std::cout << "ref: " << m1.ref() << std::endl;
  std::cout << "ref: " << m2.ref() << std::endl;
  std::cout << "ref: " << m3.ref() << std::endl;

  vec1d<float> v1 = m2.row(3);
  
  std::cout << "v1 : " << v1;
  std::cout << "ref: " << v1.ref() << std::endl;

  vec1d<float> v2 = m2.col(3);
  
  std::cout << "v2 : " << v2;
  std::cout << "ref: " << v2.ref() << std::endl;
   
  std::cout << "===============" << std::endl;
  m2 = m1.rowRange(3,5);
  
  std::cout << "m1 : " << m1;
  std::cout << "ref: " << m1.ref() << std::endl;
  std::cout << "m2 : " << m2;
  std::cout << "ref: " << m2.ref() << std::endl;

  m2 = m1.colRange(3,6);
  
  std::cout << "m1 : " << m1;
  std::cout << "ref: " << m1.ref() << std::endl;
  std::cout << "m2 : " << m2;
  std::cout << "ref: " << m2.ref() << std::endl;

  m2 = m1.subMat(0, 0, 2, 3);

  std::cout << "m2 : " << m2;
  std::cout << "ref: " << m2.ref() << std::endl;
    
}

int main()
{
  testMat2d();

  return 0;
}
