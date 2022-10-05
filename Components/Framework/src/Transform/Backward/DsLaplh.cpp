/**
 * @file DsLaplh.cpp
 * @brief Source of the backward transform operator Backard::DsLaplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/DsLaplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string DsLaplh::sTag()
   {
      return "Bwd::DsLaplh";
   }

   std::string DsLaplh::sFormatted()
   {
      return "Backard::DsLaplh";
   }

   DsLaplh::DsLaplh()
      : IRegisterId<DsLaplh>(DsLaplh::sTag(), DsLaplh::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
