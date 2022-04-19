/**
 * @file Delta.cpp
 * @brief Source of the Delta nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Delta.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Delta::sTag()
   {
      return "delta";
   }

   std::string Delta::sFormatted()
   {
      return "Delta";
   }

   Delta::Delta(const MHDFloat value)
      : IRegisterId<Delta>(value, Delta::sTag(), Delta::sFormatted())
   {
   }

}
}
