#include "../src/vec1d.hpp"

int main()
{
  // create v1 with null data
  vec1d<int> v1(3);
  v1.print();

  // create v2 and use 4 to init
  vec1d<int> v2(3, 4);
  v2.print();

  // assigment
  v1 = v2;
  v1.print();

  v1[2] = 5;
  v1.print();
  v2.print();

  // init by using array
  int array[] = {2,3,4,5,65};
  vec1d<int> v3(5, array);
  v3.print();

  // slice
  v1 = v3(1,4);
  v1[3] = 222;
  v1.print();
  v3.print();

  v1 = v3.range(1,4);
  v1.print();

  return 0;
}
