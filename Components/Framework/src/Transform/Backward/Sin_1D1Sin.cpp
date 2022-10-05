/**
 * @file Sin_1D1Sin.cpp
 * @brief Source of the backward transform operator Backard::Sin_1D1Sin
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1D1Sin.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1D1Sin::sTag()
   {
      return "Bwd::Sin_1D1Sin";
   }

   std::string Sin_1D1Sin::sFormatted()
   {
      return "Backard::Sin_1D1Sin";
   }

   Sin_1D1Sin::Sin_1D1Sin()
      : IRegisterId<Sin_1D1Sin>(Sin_1D1Sin::sTag(), Sin_1D1Sin::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
