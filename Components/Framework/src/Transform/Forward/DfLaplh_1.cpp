/**
 * @file DfLaplh_1.cpp
 * @brief Source of the forward transform operator Forward::DfLaplh_1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/DfLaplh_1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string DfLaplh_1::sTag()
   {
      return "Fwd::DfLaplh_1";
   }

   std::string DfLaplh_1::sFormatted()
   {
      return "Forward::DfLaplh_1";
   }

   DfLaplh_1::DfLaplh_1()
      : IRegisterId<DfLaplh_1>(DfLaplh_1::sTag(), DfLaplh_1::sFormatted())
   {
   }

} // Forward
} // Transform
} // QuICC
