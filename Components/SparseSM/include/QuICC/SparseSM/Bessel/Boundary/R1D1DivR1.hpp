/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of the boundary value for toroidal stress-free for Bessel polynomials
 */

#ifndef QUICC_SPARSESM_BESSEL_BOUNDARY_R1D1DIVR1_HPP
#define QUICC_SPARSESM_BESSEL_BOUNDARY_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/// Namespace for Boundary tau lines for Bessel polynomials
namespace Boundary {

   /**
    * @brief Implementation of the boundary value for toroidal stress-free for Bessel polynomial
    */
   class R1D1DivR1: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param l       Harmonic degree l
          */
         R1D1DivR1(const BesselKind type, const int l);

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
} // Bessel
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_BESSEL_BOUNDARY_R1D1DIVR1_HPP
