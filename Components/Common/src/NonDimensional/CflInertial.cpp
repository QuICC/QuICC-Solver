/**
 * @file CflInertial.cpp
 * @brief Source of the Inertial CFL nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/CflInertial.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string CflInertial::sTag()
   {
      return "cfl_inertial";
   }

   std::string CflInertial::sFormatted()
   {
      return "Inertial CFL";
   }

   CflInertial::CflInertial(const MHDFloat value)
      : IRegisterId<CflInertial>(value, CflInertial::sTag(), CflInertial::sFormatted())
   {
   }

}
}
