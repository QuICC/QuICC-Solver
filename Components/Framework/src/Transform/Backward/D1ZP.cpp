/**
 * @file D1ZP.cpp
 * @brief Source of the backward transform operator Backard::D1ZP
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1ZP.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1ZP::sTag()
   {
      return "Bwd::D1ZP";
   }

   std::string D1ZP::sFormatted()
   {
      return "Backard::D1ZP";
   }

   D1ZP::D1ZP()
      : IRegisterId<D1ZP>(D1ZP::sTag(), D1ZP::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
