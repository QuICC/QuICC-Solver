/** 
 * @file QuICCTimer.cpp
 * @brief Source of static execution timer
 */

// System includes
//

// Project includes
//
#include "QuICC/QuICCTimer.hpp"

namespace QuICC {

   ExecutionTimer& QuICCTimer()
   {
      static ExecutionTimer timer;

      return timer;
   }

}
