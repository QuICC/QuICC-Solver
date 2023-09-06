/**
 * @file InsulatingShell.hpp
 * @brief Implementation of the boundary insulating condition in a shell for Chebyshev polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGSHELL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGSHELL_HPP

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
    * @brief Implementation of the boundary insulating condition in a shell for Chebyshev polynomial
    */
   class InsulatingShell: public ICondition
   {
      public:
         /**
          * @brief Constructor for given position
          *
          * @param lower Lower bound of y
          * @param upper Upper bound of y
          * @param pos   Position of the boundary
          * @param l     Harmonic degree
          */
         InsulatingShell(const Scalar_t lower, const Scalar_t upper, const Position position, const int l);

         /**
          * @brief Destructor
          */
         ~InsulatingShell() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN       Highest polynomial
          */
         ACoeff_t compute(const int maxN);

      private:
         /**
          * @brief Harmonic degree
          */
         int mL;
   };

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGSHELL_HPP
