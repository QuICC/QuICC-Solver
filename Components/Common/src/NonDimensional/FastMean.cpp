/**
 * @file FastMean.cpp
 * @brief Source of the fast mean nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/FastMean.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string FastMean::sTag()
   {
      return "fast_mean";
   }

   std::string FastMean::sFormatted()
   {
      return "fast mean";
   }

   FastMean::FastMean(const MHDFloat value)
      : IRegisterId<FastMean>(value, FastMean::sTag(), FastMean::sFormatted())
   {
   }

}
}
