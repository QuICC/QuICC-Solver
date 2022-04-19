/**
 * @file ExplicitNextstep.cpp
 * @brief Source of the ExplicitNextstep ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string ExplicitNextstep::sTag()
   {
      return "explicit_nextstep";
   }

   std::string ExplicitNextstep::sFormatted()
   {
      return "ExplicitNextstep";
   }

   ExplicitNextstep::ExplicitNextstep()
      : IRegisterId<ExplicitNextstep>(ExplicitNextstep::sTag(), ExplicitNextstep::sFormatted())
   {
   }

}
}
