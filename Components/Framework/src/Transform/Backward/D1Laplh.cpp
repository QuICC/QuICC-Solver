/**
 * @file D1Laplh.cpp
 * @brief Source of the backward transform operator Backard::D1Laplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1Laplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1Laplh::sTag()
   {
      return "Bwd::D1Laplh";
   }

   std::string D1Laplh::sFormatted()
   {
      return "Backard::D1Laplh";
   }

   D1Laplh::D1Laplh()
      : IRegisterId<D1Laplh>(D1Laplh::sTag(), D1Laplh::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
