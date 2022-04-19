/** 
 * @file TimerTools.cpp
 * @brief Source of the timer tools
 */

// System includes
//

// External includes
//

// Class include
//
#include "TimerTools.hpp"

// Project includes
//

namespace QuICC {

   double TimerTools::reset(ITimer& rTimer)
   {
      // Stop the timer
      rTimer.stop();

      // Get elapsed time
      double tmp = rTimer.time();

      // Restart timer
      rTimer.start();

      // return elapsed time
      return tmp;
   }

}
