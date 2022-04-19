/**
 * @file RRatio.cpp
 * @brief Source of the R ratio nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/RRatio.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string RRatio::sTag()
   {
      return "rratio";
   }

   std::string RRatio::sFormatted()
   {
      return "R ratio";
   }

   RRatio::RRatio(const MHDFloat value)
      : IRegisterId<RRatio>(value, RRatio::sTag(), RRatio::sFormatted())
   {
   }

}
}
