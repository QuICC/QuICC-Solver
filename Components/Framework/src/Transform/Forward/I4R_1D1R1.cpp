/**
 * @file I4R_1D1R1.cpp
 * @brief Source of the forward transform operator Forward::I4R_1D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4R_1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4R_1D1R1::sTag()
   {
      return "Fwd::I4R_1D1R1";
   }

   std::string I4R_1D1R1::sFormatted()
   {
      return "Forward::I4R_1D1R1";
   }

   I4R_1D1R1::I4R_1D1R1()
      : IRegisterId<I4R_1D1R1>(I4R_1D1R1::sTag(), I4R_1D1R1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
