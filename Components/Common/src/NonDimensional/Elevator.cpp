/**
 * @file Elevator.cpp
 * @brief Source of the Elevator nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Elevator.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Elevator::sTag()
   {
      return "elevator";
   }

   std::string Elevator::sFormatted()
   {
      return "Elevator";
   }

   Elevator::Elevator(const MHDFloat value)
      : IRegisterId<Elevator>(value, Elevator::sTag(), Elevator::sFormatted())
   {
   }

}
}
