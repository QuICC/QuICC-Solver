/**
 * @file Lambda.cpp
 * @brief Source of the Lambda nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Lambda.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Lambda::sTag()
   {
      return "lambda";
   }

   std::string Lambda::sFormatted()
   {
      return "Lambda";
   }

   Lambda::Lambda(const MHDFloat value)
      : IRegisterId<Lambda>(value, Lambda::sTag(), Lambda::sFormatted())
   {
   }

}
}
