/**
 * @file D3.cpp
 * @brief Source of the backward transform operator Backard::D3
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D3.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D3::sTag()
   {
      return "Bwd::D3";
   }

   std::string D3::sFormatted()
   {
      return "Backard::D3";
   }

   D3::D3()
      : IRegisterId<D3>(D3::sTag(), D3::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
