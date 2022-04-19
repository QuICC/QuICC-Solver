/**
 * @file FieldToRhs.cpp
 * @brief Source of the FieldToRhs ModelOperatorBoundary
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperatorBoundary {

   std::string FieldToRhs::sTag()
   {
      return "field_to_rhs";
   }

   std::string FieldToRhs::sFormatted()
   {
      return "FieldToRhs";
   }

   FieldToRhs::FieldToRhs()
      : IRegisterId<FieldToRhs>(FieldToRhs::sTag(), FieldToRhs::sFormatted())
   {
   }

}
}
