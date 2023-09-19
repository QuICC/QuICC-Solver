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
#include "Types/Precision.hpp"
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
         JacobiRule(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Compute the quadrature
          */
         void computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size);

      protected:
         /**
          * @brief Get p polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   p(const internal::MHDLong xi, const int diff);

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   q(const internal::MHDLong xi, const int diff);

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   r(const int n, const int diff);

         /**
          * @brief Compute polynomial value
          */
         virtual internal::MHDLong   u(const internal::MHDLong x, const int n);

         /**
          * @brief Compute first derivative value
          */
         virtual internal::MHDLong  du(const internal::MHDLong x, const int n);

         /**
          * @brief Estimate first zero
          *
          * Theorem 1 from D. Dimitrov, G. Nikolov, Sharp bounds for the extreme zeros of classical orthogonal polynomials, Journal of Approximation Theory, 2010
          */
         internal::MHDLong  estimateFirstZero(const int n);

      private:
         /**
          * @brief Alpha
          */
         internal::MHDFloat mAlpha;

         /**
          * @brief Beta
          */
         internal::MHDFloat mBeta;
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_JACOBIRULE_HPP
