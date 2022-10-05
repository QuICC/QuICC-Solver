/**
 * @file I4R_1D1R1ZI2.cpp
 * @brief Source of the forward transform operator Forward::I4R_1D1R1ZI2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4R_1D1R1ZI2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4R_1D1R1ZI2::sTag()
   {
      return "Fwd::I4R_1D1R1ZI2";
   }

   std::string I4R_1D1R1ZI2::sFormatted()
   {
      return "Forward::I4R_1D1R1ZI2";
   }

   I4R_1D1R1ZI2::I4R_1D1R1ZI2()
      : IRegisterId<I4R_1D1R1ZI2>(I4R_1D1R1ZI2::sTag(), I4R_1D1R1ZI2::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
