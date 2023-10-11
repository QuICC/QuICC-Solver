/**
 * @file LegendreRule.hpp
 * @brief Implementation of a Legendre quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_LEGENDRERULE_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_LEGENDRERULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Quadrature/PrueferAlgorithm.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of a Legendre quadrature rule
    */
   class LegendreRule: public PrueferAlgorithm
   {
      public:
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
         virtual Internal::MHDLong   p(const Internal::MHDLong xi, const int diff) override;

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   q(const Internal::MHDLong xi, const int diff) override;

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   r(const int size, const int diff) override;

         /**
          * @brief Node estimate
          *
          * @param k Index of the node to estimate
          */
         Internal::MHDLong   estimateNode(const int k, const int size);

         /**
          * @brief Compute polynomial value at 0
          */
         Internal::MHDLong   zeroPoly(const int size);

         /**
          * @brief Compute first derivative value at 0
          */
         Internal::MHDLong   zeroDiff(const int size);

      private:
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_LEGENDRERULE_HPP
