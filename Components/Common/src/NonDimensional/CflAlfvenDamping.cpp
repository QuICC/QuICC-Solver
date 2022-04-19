/**
 * @file CflAlfvenDamping.cpp
 * @brief Source of the Alfven damping for CFL nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string CflAlfvenDamping::sTag()
   {
      return "cfl_alfven_damping";
   }

   std::string CflAlfvenDamping::sFormatted()
   {
      return "Alfven damping for CFL";
   }

   CflAlfvenDamping::CflAlfvenDamping(const MHDFloat value)
      : IRegisterId<CflAlfvenDamping>(value, CflAlfvenDamping::sTag(), CflAlfvenDamping::sFormatted())
   {
   }

}
}
