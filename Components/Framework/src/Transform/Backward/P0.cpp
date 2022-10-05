/**
 * @file P0.cpp
 * @brief Source of the backward transform operator Backard::P0
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/P0.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string P0::sTag()
   {
      return "Bwd::P0";
   }

   std::string P0::sFormatted()
   {
      return "Backard::P0";
   }

   P0::P0()
      : IRegisterId<P0>(P0::sTag(), P0::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
