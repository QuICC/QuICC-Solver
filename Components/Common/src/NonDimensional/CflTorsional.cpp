/**
 * @file CflTorsional.cpp
 * @brief Source of the Torsional CFL nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/CflTorsional.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string CflTorsional::sTag()
   {
      return "cfl_torsional";
   }

   std::string CflTorsional::sFormatted()
   {
      return "Torsional CFL";
   }

   CflTorsional::CflTorsional(const MHDFloat value)
      : IRegisterId<CflTorsional>(value, CflTorsional::sTag(), CflTorsional::sFormatted())
   {
   }

}
}
