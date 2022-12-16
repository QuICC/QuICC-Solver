/** 
 * @file Tools.cpp
 * @brief Source of the tools for full sphere Worland operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Tools.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   WorlandKind Tools::identifyBasis(const Scalar_t alpha, const Scalar_t dBeta)
   {
      WorlandKind type;

      if(alpha == MHD_MP(-0.5))
      {
         type = WorlandKind::CHEBYSHEV;
      } else if(alpha == MHD_MP(0.0))
      {
         if(dBeta == MHD_MP(-0.5))
         {
            type = WorlandKind::LEGENDRE;
         } else if(dBeta == MHD_MP(0.0))
         {
            type = WorlandKind::CYLENERGY;
         } else if(dBeta == MHD_MP(0.5))
         {
            type = WorlandKind::SPHENERGY;
         } else
         {
            throw std::logic_error("Unknown Legendre type Worland for SparseSM");
         }
      } else
      {
         throw std::logic_error("Unknown Worland type for SparseSM");
      }

      return type;
   }

}
}
}
