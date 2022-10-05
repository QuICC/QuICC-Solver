/**
 * @file Laplh_1Sin_1Dphi.cpp
 * @brief Source of the forward transform operator Forward::Laplh_1Sin_1Dphi
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1Sin_1Dphi.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1Sin_1Dphi::sTag()
   {
      return "Fwd::Laplh_1Sin_1Dphi";
   }

   std::string Laplh_1Sin_1Dphi::sFormatted()
   {
      return "Forward::Laplh_1Sin_1Dphi";
   }

   Laplh_1Sin_1Dphi::Laplh_1Sin_1Dphi()
      : IRegisterId<Laplh_1Sin_1Dphi>(Laplh_1Sin_1Dphi::sTag(), Laplh_1Sin_1Dphi::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
