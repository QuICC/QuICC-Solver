/**
 * @file I6R_1D1R1.cpp
 * @brief Source of the forward transform operator Forward::I6R_1D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6R_1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6R_1D1R1::sTag()
   {
      return "Fwd::I6R_1D1R1";
   }

   std::string I6R_1D1R1::sFormatted()
   {
      return "Forward::I6R_1D1R1";
   }

   I6R_1D1R1::I6R_1D1R1()
      : IRegisterId<I6R_1D1R1>(I6R_1D1R1::sTag(), I6R_1D1R1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
