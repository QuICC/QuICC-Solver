/**
 * @file D1.cpp
 * @brief Source of the backward transform operator Backward::D1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1::sTag()
   {
      return "Bwd::D1";
   }

   std::string D1::sFormatted()
   {
      return "Backward::D1";
   }

   D1::D1()
      : IRegisterId<D1>(D1::sTag(), D1::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
