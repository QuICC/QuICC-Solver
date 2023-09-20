/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of the boundary value of r D 1/r for Chebyshev linear map polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_R1D1DIVR1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   /**
    * @brief Implementation of the boundary value of r D 1/r for Worland polynomial
    */
   class R1D1DivR1: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param lower Lower bound of y
          * @param upper Upper bound of y
          * @param pos   Position of the boundary
          */
         R1D1DivR1(const Scalar_t lower, const Scalar_t upper, const Position pos);

         /**
          * @brief Destructor
          */
         ~R1D1DivR1() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN Highest polynomial
          */
         ACoeff_t compute(const int maxN);

      private:
   };

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_R1D1DIVR1_HPP
