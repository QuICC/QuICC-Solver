/**
 * @file Undefined.cpp
 * @brief Source of the Undefined physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Undefined.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Undefined::sTag()
   {
      return "undefined";
   }

   std::string Undefined::sFormatted()
   {
      return "Undefined";
   }

   Undefined::Undefined()
      : IRegisterId<Undefined>(Undefined::sTag(), Undefined::sFormatted())
   {
   }

}
}
