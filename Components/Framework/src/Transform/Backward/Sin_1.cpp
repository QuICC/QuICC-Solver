/**
 * @file Sin_1.cpp
 * @brief Source of the backward transform operator Backard::Sin_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1::sTag()
   {
      return "Bwd::Sin_1";
   }

   std::string Sin_1::sFormatted()
   {
      return "Backard::Sin_1";
   }

   Sin_1::Sin_1()
      : IRegisterId<Sin_1>(Sin_1::sTag(), Sin_1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
