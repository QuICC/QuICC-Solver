/**
 * @file DfLaplh.cpp
 * @brief Source of the backward transform operator Backard::DfLaplh
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/DfLaplh.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string DfLaplh::sTag()
   {
      return "Bwd::DfLaplh";
   }

   std::string DfLaplh::sFormatted()
   {
      return "Backard::DfLaplh";
   }

   DfLaplh::DfLaplh()
      : IRegisterId<DfLaplh>(DfLaplh::sTag(), DfLaplh::sFormatted())
   {
   }

} // Backward
} // Transform
} // QuICC
