/**
 * @file TimerTools.hpp
 * @brief Implementation of timer tools 
 */

#ifndef QUICC_TIMERTOOLS_HPP
#define QUICC_TIMERTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "ITimer.hpp"

namespace QuICC {

   /**
    * @brief Implementation of timer tools
    */
   class TimerTools
   {
      public:
         /**
          * @brief Reset timer (stop and restart) and return elapsed time
          */
         static double reset(ITimer& rTimer);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         TimerTools();

         /**
          * @brief Destructor
          */
         ~TimerTools();

   };

}

#endif // QUICC_TIMERTOOLS_HPP
