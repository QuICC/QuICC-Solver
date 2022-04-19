/**
 * @file Taylor.cpp
 * @brief Source of the Taylor nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Taylor.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Taylor::sTag()
   {
      return "taylor";
   }

   std::string Taylor::sFormatted()
   {
      return "Taylor";
   }

   Taylor::Taylor(const MHDFloat value)
      : IRegisterId<Taylor>(value, Taylor::sTag(), Taylor::sFormatted())
   {
   }

}
}
