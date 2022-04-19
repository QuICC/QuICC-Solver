/**
 * @file Boundary.cpp
 * @brief Source of the Boundary ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/Boundary.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string Boundary::sTag()
   {
      return "boundary";
   }

   std::string Boundary::sFormatted()
   {
      return "Boundary";
   }

   Boundary::Boundary()
      : IRegisterId<Boundary>(Boundary::sTag(), Boundary::sFormatted())
   {
   }

}
}
