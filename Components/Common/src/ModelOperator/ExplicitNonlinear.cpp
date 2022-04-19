/**
 * @file ExplicitNonlinear.cpp
 * @brief Source of the ExplicitNonlinear ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string ExplicitNonlinear::sTag()
   {
      return "explicit_nonlinear";
   }

   std::string ExplicitNonlinear::sFormatted()
   {
      return "ExplicitNonlinear";
   }

   ExplicitNonlinear::ExplicitNonlinear()
      : IRegisterId<ExplicitNonlinear>(ExplicitNonlinear::sTag(), ExplicitNonlinear::sFormatted())
   {
   }

}
}
