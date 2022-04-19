/**
 * @file SerialDebugger.hpp
 * @brief Implementation of a very simple serial debugger 
 */

#ifndef QUICC_DEBUG_SERIALDEBUGGER_HPP
#define QUICC_DEBUG_SERIALDEBUGGER_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "Timers/TimerMacro.h"

namespace QuICC {

namespace Debug {

   /**
    * @brief Implementation of a very simple serial debugger
    */
   class SerialDebugger
   {
      public:
         /**
          * @brief Debug message
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          */
         static void msg(const std::string& msg, const int tabs);

         /**
          * @brief Debug message when entering function
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          */
         static void enter(const std::string& msg, const int tabs);

         /**
          * @brief Debug message when leaving function
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          */
         static void leave(const std::string& msg, const int tabs);

         /**
          * @brief Debug message with numerical value
          *
          * @param msg Message to print
          * @param tabs Number of tab characters
          * @param value Value to add at end of message
          */
         static void showValue(const std::string& msg, const int tabs, const MHDFloat value);

         /**
          * @brief Timer to use for debugging
          */
         static TimerMacro timer;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         SerialDebugger();

         /**
          * @brief Destructor
          */
         ~SerialDebugger();
   };

}
}

#endif // QUICC_DEBUG_SERIALDEBUGGER_HPP
