/**
 * @file Beta.cpp
 * @brief Source of the Beta nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Beta.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Beta::sTag()
   {
      return "beta";
   }

   std::string Beta::sFormatted()
   {
      return "Beta";
   }

   Beta::Beta(const MHDFloat value)
      : IRegisterId<Beta>(value, Beta::sTag(), Beta::sFormatted())
   {
   }

}
}
