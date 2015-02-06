/*!
 * \file global.hpp
 * \breif some commonly used definitions.
 *
 * \author Shaoguang Cheng. From Xi'an, China
 * \date   1/2/2015
 */

#ifndef __GLOBAL__
#define __GLOBAL__

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#ifndef DEBUGMSG
#define DEBUGMSG(msg) std::cout << "line: " << __LINE__	\
  << ", file: " << __FILE__ \
  << ", message: " << msg << std::endl
#endif // end of debugmsg

class useCount
{
 private:
  int *n;
  
 public:
 useCount() :  n(new int(1)) {}
 useCount(const useCount& u) : n(u.n){++*n;}
  ~useCount(){
    if(--*n == 0) 
      delete n;
  }

  bool only() const {
    if(*n == 1)
      return true;
    else
      return false;
  }

  bool reattach(const useCount& u){
    ++*u.n;
    if(1 == --*n){
      delete n;
      n = NULL;
      n = u.n;
      return true;
    }

    n = u.n;
    return false;
  }

  int getCount() const{
    return *n;
  }

 private:
  useCount& operator = (const useCount& u);
};

#endif // end of global
