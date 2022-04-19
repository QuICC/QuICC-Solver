/** 
 * @file WorlandLegendreRule.cpp
 * @brief Source of the Worland Legendre quadrature
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void WorlandLegendreRule::computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);
      int jsize = 2*size;
      internal::Array jgrid(jsize);
      internal::Array jweights(jsize);
      LegendreRule::computeQuadrature(jgrid, jweights, jsize);

      igrid = jgrid.bottomRows(size);
      iweights = jweights.bottomRows(size);
   }

}
}
}
