/**
 * @file StressFreePolAnelastic.hpp
 * @brief Implementation of the boundary value of r D(* / r) - F  for Chebyshev linear map polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_STRESSFREEPOLANELASTIC_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_STRESSFREEPOLANELASTIC_HPP

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
    * @brief Implementation of the boundary value of r D(* / r) - F  for Worland polynomial
    */
   class StressFreePolAnelastic: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param lower Lower bound of y
          * @param upper Upper bound of y
          * @param pos   Position of the boundary
          * @param Fb    Boundary value of the radial field
          */
         StressFreePolAnelastic(const Scalar_t lower, const Scalar_t upper, const Position pos, const Internal::MHDFloat Fb);

         /**
          * @brief Destructor
          */
         ~StressFreePolAnelastic() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN Highest polynomial
          */
         ACoeff_t compute(const int maxN);

      private:
         /**
          * @brief Radial field boundary value
          */
         Internal::MHDFloat mFb;
   };

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_STRESSFREEPOLANELASTIC_HPP
