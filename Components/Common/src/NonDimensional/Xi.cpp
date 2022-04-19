/**
 * @file Xi.cpp
 * @brief Source of the Xi nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Xi.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Xi::sTag()
   {
      return "xi";
   }

   std::string Xi::sFormatted()
   {
      return "Xi";
   }

   Xi::Xi(const MHDFloat value)
      : IRegisterId<Xi>(value, Xi::sTag(), Xi::sFormatted())
   {
   }

}
}
