/** 
 * @file SerialDebugger.cpp
 * @brief Source of the serial debugger implementation
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "QuICC/Debug/SerialDebugger.hpp"

// Project includes
//

namespace QuICC {

namespace Debug {

   TimerMacro SerialDebugger::timer = TimerMacro();

   void SerialDebugger::msg(const std::string& msg, const int tabs)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "DEBUG " << msg << std::endl;
   }

   void SerialDebugger::enter(const std::string& msg, const int tabs)
   {
      SerialDebugger::msg("ENTERING: " + msg, tabs);
   }

   void SerialDebugger::leave(const std::string& msg, const int tabs)
   {
      SerialDebugger::msg("LEAVING: " + msg, tabs);
   }

   void SerialDebugger::showValue(const std::string& msg, const int tabs, const MHDFloat value)
   {
      for(int i = 0; i < tabs; i++)
      {
         std::cerr << "\t";
      }
      std::cerr << "DEBUG " << msg << value << std::endl;
   }

}
}
