/**
 * @file Sin_1Laplh.cpp
 * @brief Source of the backward transform operator Backard::Sin_1Laplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1Laplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1Laplh::sTag()
   {
      return "Bwd::Sin_1Laplh";
   }

   std::string Sin_1Laplh::sFormatted()
   {
      return "Backard::Sin_1Laplh";
   }

   Sin_1Laplh::Sin_1Laplh()
      : IRegisterId<Sin_1Laplh>(Sin_1Laplh::sTag(), Sin_1Laplh::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
