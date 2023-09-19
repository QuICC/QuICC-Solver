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
#include "Types/Precision.hpp"
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
         void computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size);

      protected:
         /**
          * @brief Get p polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   p(const internal::MHDLong xi, const int diff) override;

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   q(const internal::MHDLong xi, const int diff) override;

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   r(const int size, const int diff) override;

         /**
          * @brief Node estimate
          *
          * @param k Index of the node to estimate
          */
         internal::MHDLong   estimateNode(const int k, const int size);

         /**
          * @brief Compute polynomial value at 0
          */
         internal::MHDLong   zeroPoly(const int size);

         /**
          * @brief Compute first derivative value at 0
          */
         internal::MHDLong   zeroDiff(const int size);

      private:
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_LEGENDRERULE_HPP
