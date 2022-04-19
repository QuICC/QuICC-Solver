/**
 * @file Prognostic.cpp
 * @brief Source of the Prognostic SolveTiming
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SolveTiming/Prognostic.hpp"

// Project includes
//

namespace QuICC {

namespace SolveTiming {

   std::string Prognostic::sTag()
   {
      return "prognostic";
   }

   std::string Prognostic::sFormatted()
   {
      return "Prognostic";
   }

   Prognostic::Prognostic()
      : IRegisterId<Prognostic>(Prognostic::sTag(), Prognostic::sFormatted())
   {
   }

}
}
