/**
 * @file ILinearMapOperator.hpp
 * @brief Implementation of the generic interface to the spherical shell with linear map dense operator
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ILINEARMAPOPERATOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ILINEARMAPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the generic interface to the spherical shell with linear map dense operator
    */
   class ILinearMapOperator: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound
          * @param upper   Upper bound
          */
         ILinearMapOperator(const int rows, const int cols, const Scalar_t  lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~ILinearMapOperator() = default;

      protected:
         /**
          * @brief Compute quadrature grid and weights
          */
         void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size) const;

         /**
          * @brief Lower bound of domain
          */
         const Scalar_t mcLower;

         /**
          * @brief  Upper bound of domain
          */
         const Scalar_t mcUpper;

      private:
   };

}
}
}
}

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_ILINEARMAPOPERATOR_HPP
