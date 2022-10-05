/**
 * @file Laplh.cpp
 * @brief Source of the backward transform operator Backard::Laplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Laplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Laplh::sTag()
   {
      return "Bwd::Laplh";
   }

   std::string Laplh::sFormatted()
   {
      return "Backard::Laplh";
   }

   Laplh::Laplh()
      : IRegisterId<Laplh>(Laplh::sTag(), Laplh::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
