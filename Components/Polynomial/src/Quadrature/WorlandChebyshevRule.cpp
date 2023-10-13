/**
 * @file WorlandChebyshevRule.cpp
 * @brief Source of the Worland Chebyshev quadrature
 */

// System includes
//

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandChebyshevRule.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void WorlandChebyshevRule::computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);

      for(int k = 0; k < size/2; k++)
      {
         Internal::MHDFloat theta = Internal::Math::PI_long*(Internal::MHDFloat(2*k + 1))/Internal::MHDFloat(4*size);
         igrid(k) = Internal::Math::sin(theta);
      }

      for(int k = size/2; k < size; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         Internal::MHDFloat theta = Internal::Math::PI_long*(Internal::MHDFloat(2*k_ - 1))/Internal::MHDFloat(4*size);
         igrid(k) = Internal::Math::cos(theta);
      }

      iweights.setConstant(Internal::Math::PI/(MHD_MP(2.0)*Internal::MHDFloat(size)));
   }

   void WorlandChebyshevRule::computeXQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);

      for(int k = 0; k < size/3; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         Internal::MHDFloat theta = Internal::Math::PI_long*(Internal::MHDFloat(2*k_ - 1))/Internal::MHDFloat(2*size);
         igrid(k) = Internal::Math::cos(theta);
      }

      for(int k = size/3; k < 2*size/3; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         Internal::MHDFloat theta = Internal::Math::PI_long*(Internal::MHDFloat(size-2*k_ + 1))/Internal::MHDFloat(2*size);
         igrid(k) = Internal::Math::sin(theta);
      }

      for(int k = 2*size/3; k < size; k++)
      {
         // Reverse grid r = (0, 1)
         int k_ = size-k;

         Internal::MHDFloat theta = Internal::Math::PI_long*(Internal::MHDFloat(2*k_ - 1))/Internal::MHDFloat(2*size);
         igrid(k) = Internal::Math::cos(theta);
      }

      iweights.setConstant(Internal::Math::PI/(MHD_MP(2.0)*Internal::MHDFloat(size)));
   }

}
}
}
