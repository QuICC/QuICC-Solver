/**
 * @file AnnulusExactStateIds.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cylindrical annulus geometries
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/AnnulusExactStateIds.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Equations {

   const MHDFloat AnnulusExactStateIds::PCOS = 99999;

   const MHDFloat AnnulusExactStateIds::PSIN = -99999;

   MHDFloat AnnulusExactStateIds::cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::cos(mode*theta);
   }

   MHDFloat AnnulusExactStateIds::sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::sin(mode*theta);
   }

   MHDFloat AnnulusExactStateIds::poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      MHDFloat val;

      if(mode == AnnulusExactStateIds::PCOS)
      {
         val = AnnulusExactStateIds::cos(amplitude,mode,Math::PI*(x-1)/2.0);
      } else if(mode == AnnulusExactStateIds::PSIN)
      {
         val = AnnulusExactStateIds::sin(amplitude,mode,Math::PI*(x-1)/2.0);
      } else
      {
         val = amplitude*std::pow(x,mode);
      }

      return val;
   }

}
}
