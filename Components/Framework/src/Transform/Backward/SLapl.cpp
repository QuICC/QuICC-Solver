/**
 * @file SLapl.cpp
 * @brief Source of the backward transform operator Backard::SLapl
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/SLapl.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string SLapl::sTag()
   {
      return "Bwd::SLapl";
   }

   std::string SLapl::sFormatted()
   {
      return "Backard::SLapl";
   }

   SLapl::SLapl()
      : IRegisterId<SLapl>(SLapl::sTag(), SLapl::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
