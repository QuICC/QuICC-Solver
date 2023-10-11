/**
 * @file D1.hpp
 * @brief Implementation of the boundary first derivative for Chebyshev polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_D1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_D1_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   /**
    * @brief Implementation of the boundary first derivative for Chebyshev polynomial
    */
   class D1: public ICondition
   {
      public:
         /**
          * @brief Constructor for given position
          *
          * @param lower Lower bound of y
          * @param upper Upper bound of y
          * @param pos   Position of the boundary
          */
         D1(const Scalar_t lower, const Scalar_t upper, const Position position);

         /**
          * @brief Destructor
          */
         ~D1() = default;

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

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_D1_HPP
