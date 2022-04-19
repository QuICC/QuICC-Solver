/**
 * @file Temperature.cpp
 * @brief Source of the Temperature physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Temperature.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Temperature::sTag()
   {
      return "temperature";
   }

   std::string Temperature::sFormatted()
   {
      return "Temperature";
   }

   Temperature::Temperature()
      : IRegisterId<Temperature>(Temperature::sTag(), Temperature::sFormatted())
   {
   }

}
}
