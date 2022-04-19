/**
 * @file CflAlfvenScale.cpp
 * @brief Source of the Alfven scale for CFL nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string CflAlfvenScale::sTag()
   {
      return "cfl_alfven_scale";
   }

   std::string CflAlfvenScale::sFormatted()
   {
      return "Alfven scale for CFL";
   }

   CflAlfvenScale::CflAlfvenScale(const MHDFloat value)
      : IRegisterId<CflAlfvenScale>(value, CflAlfvenScale::sTag(), CflAlfvenScale::sFormatted())
   {
   }

}
}
