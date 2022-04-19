/**
 * @file Chi.cpp
 * @brief Source of the Chi nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Chi.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Chi::sTag()
   {
      return "chi";
   }

   std::string Chi::sFormatted()
   {
      return "Chi";
   }

   Chi::Chi(const MHDFloat value)
      : IRegisterId<Chi>(value, Chi::sTag(), Chi::sFormatted())
   {
   }

}
}
