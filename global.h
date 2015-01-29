#ifndef __GLOBAL__
#define __GLOBAL__

#include <iostream>

#ifndef DEBUGMSG
#define DEBUGMSG std::cout << "line: " << __LINE__ \
  << ", file: " << __FILE__ \
  << ", message: " << msg << endl
#endif // end of debugmsg


#endif // end of global
