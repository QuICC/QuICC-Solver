/**
 * @file CylinderExactStateIds.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a cylinder
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/CylinderExactStateIds.hpp"

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Equations {

   const MHDFloat CylinderExactStateIds::PCOS = 99999;

   const MHDFloat CylinderExactStateIds::PSIN = -99999;

   MHDFloat CylinderExactStateIds::cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::cos(mode*theta);
   }

   MHDFloat CylinderExactStateIds::sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::sin(mode*theta);
   }

   MHDFloat CylinderExactStateIds::poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      MHDFloat val;

      if(mode == CylinderExactStateIds::PCOS)
      {
         val = CylinderExactStateIds::cos(amplitude,mode,Math::PI*(x-1)/2.0);
      } else if(mode == CylinderExactStateIds::PSIN)
      {
         val = CylinderExactStateIds::sin(amplitude,mode,Math::PI*(x-1)/2.0);
      } else
      {
         val = amplitude*std::pow(x,mode);
      }

      return val;
   }

   MHDFloat CylinderExactStateIds::chebyshev(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      return amplitude*std::cos(mode*std::acos(x));
   }

}
}
