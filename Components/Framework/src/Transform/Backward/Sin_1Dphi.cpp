/**
 * @file Sin_1Dphi.cpp
 * @brief Source of the backward transform operator Backard::Sin_1Dphi
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1Dphi.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1Dphi::sTag()
   {
      return "Bwd::Sin_1Dphi";
   }

   std::string Sin_1Dphi::sFormatted()
   {
      return "Backard::Sin_1Dphi";
   }

   Sin_1Dphi::Sin_1Dphi()
      : IRegisterId<Sin_1Dphi>(Sin_1Dphi::sTag(), Sin_1Dphi::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
