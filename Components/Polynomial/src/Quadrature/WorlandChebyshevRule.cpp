/** 
 * @file WorlandChebyshevRule.cpp
 * @brief Source of the Worland Chebyshev quadrature
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Quadrature/WorlandChebyshevRule.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void WorlandChebyshevRule::computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);

      for(int k = 0; k < size/2; k++)
      {
         internal::MHDFloat theta = Precision::PI_long*(internal::MHDFloat(2*k + 1))/internal::MHDFloat(4*size);
         igrid(k) = precision::sin(theta);
      }

      for(int k = size/2; k < size; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         internal::MHDFloat theta = Precision::PI_long*(internal::MHDFloat(2*k_ - 1))/internal::MHDFloat(4*size);
         igrid(k) = precision::cos(theta);
      }

      iweights.setConstant(Precision::PI/(MHD_MP(2.0)*internal::MHDFloat(size)));
   }

   void WorlandChebyshevRule::computeXQuadrature(internal::Array& igrid, internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);

      for(int k = 0; k < size/3; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         internal::MHDFloat theta = Precision::PI_long*(internal::MHDFloat(2*k_ - 1))/internal::MHDFloat(2*size);
         igrid(k) = precision::cos(theta);
      }

      for(int k = size/3; k < 2*size/3; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         internal::MHDFloat theta = Precision::PI_long*(internal::MHDFloat(size-2*k_ + 1))/internal::MHDFloat(2*size);
         igrid(k) = precision::sin(theta);
      }

      for(int k = 2*size/3; k < size; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         internal::MHDFloat theta = Precision::PI_long*(internal::MHDFloat(2*k_ - 1))/internal::MHDFloat(2*size);
         igrid(k) = precision::cos(theta);
      }

      iweights.setConstant(Precision::PI/(MHD_MP(2.0)*internal::MHDFloat(size)));
   }

}
}
}
