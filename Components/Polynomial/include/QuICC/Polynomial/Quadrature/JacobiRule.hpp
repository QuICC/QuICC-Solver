/**
 * @file JacobiRule.hpp
 * @brief Implementation of a Jacobi quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_JACOBIRULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_JACOBIRULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/Quadrature/PrueferAlgorithm.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of a Jacobi quadrature rule
    */
   class JacobiRule: public PrueferAlgorithm
   {
      public:
         /**
          * @brief Empty constructor
          */
         JacobiRule(const Internal::MHDFloat alpha, const Internal::MHDFloat beta);

         /**
          * @brief Compute the quadrature
          */
         void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size);

      protected:
         /**
          * @brief Get p polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   p(const Internal::MHDLong xi, const int diff);

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   q(const Internal::MHDLong xi, const int diff);

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   r(const int n, const int diff);

         /**
          * @brief Compute polynomial value
          */
         virtual Internal::MHDLong   u(const Internal::MHDLong x, const int n);

         /**
          * @brief Compute first derivative value
          */
         virtual Internal::MHDLong  du(const Internal::MHDLong x, const int n);

         /**
          * @brief Estimate first zero
          *
          * Theorem 1 from D. Dimitrov, G. Nikolov, Sharp bounds for the extreme zeros of classical orthogonal polynomials, Journal of Approximation Theory, 2010
          */
         Internal::MHDLong  estimateFirstZero(const int n);

      private:
         /**
          * @brief Alpha
          */
         Internal::MHDFloat mAlpha;

         /**
          * @brief Beta
          */
         Internal::MHDFloat mBeta;
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_JACOBIRULE_HPP
