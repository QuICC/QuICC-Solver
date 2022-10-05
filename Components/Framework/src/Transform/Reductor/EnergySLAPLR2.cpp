/**
 * @file EnergySLAPLR2.cpp
 * @brief Source of the reductor transform operator Reductor::EnergySLAPLR2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergySLAPLR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergySLAPLR2::sTag()
   {
      return "Red::EnergySLAPLR2";
   }

   std::string EnergySLAPLR2::sFormatted()
   {
      return "Reductor::EnergySLAPLR2";
   }

   EnergySLAPLR2::EnergySLAPLR2()
      : IRegisterId<EnergySLAPLR2>(EnergySLAPLR2::sTag(), EnergySLAPLR2::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
