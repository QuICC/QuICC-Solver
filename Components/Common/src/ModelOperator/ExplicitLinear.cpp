/**
 * @file ExplicitLinear.cpp
 * @brief Source of the ExplicitLinear ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/ExplicitLinear.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string ExplicitLinear::sTag()
   {
      return "explicit_linear";
   }

   std::string ExplicitLinear::sFormatted()
   {
      return "ExplicitLinear";
   }

   ExplicitLinear::ExplicitLinear()
      : IRegisterId<ExplicitLinear>(ExplicitLinear::sTag(), ExplicitLinear::sFormatted())
   {
   }

}
}
