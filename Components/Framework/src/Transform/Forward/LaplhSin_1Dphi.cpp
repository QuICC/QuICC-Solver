/**
 * @file LaplhSin_1Dphi.cpp
 * @brief Source of the forward transform operator Forward::LaplhSin_1Dphi
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhSin_1Dphi.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhSin_1Dphi::sTag()
   {
      return "Fwd::LaplhSin_1Dphi";
   }

   std::string LaplhSin_1Dphi::sFormatted()
   {
      return "Forward::LaplhSin_1Dphi";
   }

   LaplhSin_1Dphi::LaplhSin_1Dphi()
      : IRegisterId<LaplhSin_1Dphi>(LaplhSin_1Dphi::sTag(), LaplhSin_1Dphi::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
