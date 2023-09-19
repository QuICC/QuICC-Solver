/**
 * @file IPlaneOperator.hpp
 * @brief Base implementation of a plane layer operator, with y = ax + b and k1, k2 Fourier modes
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_IPLANEOPERATOR_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_IPLANEOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Base implementation of a plane layer operator, with y = ax + b and k1, k2 Fourier modes
    */
   class IPlaneOperator: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of cols
          * @param lower   Lower bound
          * @param upper   Upper bound
          * @param k1      First Fourier mode
          * @param k2      Second Fourier mode
          */
         IPlaneOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2);

         /**
          * @brief Destructor
          */
         virtual ~IPlaneOperator() = default;

      protected:
         /**
          * @brief Fourier mode in first direction
          */
         Scalar_t k1() const;

         /**
          * @brief Fourier mode in second direction
          */
         Scalar_t k2() const;

      private:
         /**
          * @brief Fourier mode in first direction
          */
         Scalar_t mK1;

         /**
          * @brief Fourier mode in second direction
          */
         Scalar_t mK2;
   };

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_IPLANEOPERATOR_HPP
