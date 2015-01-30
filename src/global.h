#ifndef __GLOBAL__
#define __GLOBAL__

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <assert.h>

#ifndef DEBUGMSG
#define DEBUGMSG(msg) std::cout << "line: " << __LINE__	\
  << ", file: " << __FILE__ \
  << ", message: " << msg << std::endl
#endif // end of debugmsg


#endif // end of global
