/**
 * @file Value.hpp
 * @brief Implementation of the boundary value for Chebyshev polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_VALUE_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_VALUE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   /**
    * @brief Implementation of the boundary value for Chebyshev polynomial
    */
   class Value: public ICondition
   {
      public:
         /**
          * @brief Constructor for given position
          *
          * @param lower Lower bound of y
          * @param upper Upper bound of y
          * @param pos Position of the boundary
          */
         Value(const Scalar_t lower, const Scalar_t upper, const Position pos);

         /**
          * @brief Destructor
          */
         ~Value() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN       Highest polynomial
          */
         ACoeff_t compute(const int maxN);

      private:
   };

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_VALUE_HPP
