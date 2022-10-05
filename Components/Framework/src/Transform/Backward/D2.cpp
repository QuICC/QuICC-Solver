/**
 * @file D2.cpp
 * @brief Source of the backward transform operator Backard::D2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D2::sTag()
   {
      return "Bwd::D2";
   }

   std::string D2::sFormatted()
   {
      return "Backard::D2";
   }

   D2::D2()
      : IRegisterId<D2>(D2::sTag(), D2::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
