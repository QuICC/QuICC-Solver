/**
 * @file P.cpp
 * @brief Source of the backward transform operator Backard::P
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string P::sTag()
   {
      return "Bwd::P";
   }

   std::string P::sFormatted()
   {
      return "Backard::P";
   }

   P::P()
      : IRegisterId<P>(P::sTag(), P::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
