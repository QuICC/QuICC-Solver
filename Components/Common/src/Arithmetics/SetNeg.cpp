/**
 * @file SetNeg.cpp
 * @brief Source of the SetNeg Arithmetics
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Arithmetics/SetNeg.hpp"

// Project includes
//

namespace QuICC {

namespace Arithmetics {

   std::string SetNeg::sTag()
   {
      return "setneg";
   }

   std::string SetNeg::sFormatted()
   {
      return "SetNeg";
   }

   SetNeg::SetNeg()
      : IRegisterId<SetNeg>(SetNeg::sTag(), SetNeg::sFormatted())
   {
   }

}
}
