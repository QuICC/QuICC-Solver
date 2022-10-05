/**
 * @file SRadLapl.cpp
 * @brief Source of the backward transform operator Backard::SRadLapl
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/SRadLapl.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string SRadLapl::sTag()
   {
      return "Bwd::SRadLapl";
   }

   std::string SRadLapl::sFormatted()
   {
      return "Backard::SRadLapl";
   }

   SRadLapl::SRadLapl()
      : IRegisterId<SRadLapl>(SRadLapl::sTag(), SRadLapl::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
