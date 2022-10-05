/**
 * @file D1R1.cpp
 * @brief Source of the backward transform operator Backard::D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1R1::sTag()
   {
      return "Bwd::D1R1";
   }

   std::string D1R1::sFormatted()
   {
      return "Backard::D1R1";
   }

   D1R1::D1R1()
      : IRegisterId<D1R1>(D1R1::sTag(), D1R1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
