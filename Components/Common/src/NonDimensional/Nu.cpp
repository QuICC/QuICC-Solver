/**
 * @file Nu.cpp
 * @brief Source of the Nu nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Nu.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Nu::sTag()
   {
      return "nu";
   }

   std::string Nu::sFormatted()
   {
      return "Nu";
   }

   Nu::Nu(const MHDFloat value)
      : IRegisterId<Nu>(value, Nu::sTag(), Nu::sFormatted())
   {
   }

}
}
