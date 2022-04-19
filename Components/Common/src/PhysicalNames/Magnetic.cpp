/**
 * @file Magnetic.cpp
 * @brief Source of the Magnetic physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Magnetic.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Magnetic::sTag()
   {
      return "magnetic";
   }

   std::string Magnetic::sFormatted()
   {
      return "Magnetic";
   }

   Magnetic::Magnetic()
      : IRegisterId<Magnetic>(Magnetic::sTag(), Magnetic::sFormatted())
   {
   }

}
}
