/**
 * @file InsulatingTopPerfcondBottomPolShell.hpp
 * @brief Implementation of the boundary mixed condition insulating top / perfectly conducting bottom in a shell for
 * Chebyshev polynomials
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGTOPPERFCONDBOTTOMPOLSHELL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGTOPPERFCONDBOTTOMPOLSHELL_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

/**
 * @brief Implementation of the boundary insulating condition in a shell for
 * Chebyshev polynomial
 */
class InsulatingTopPerfcondBottomPolShell : public ICondition
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
   InsulatingTopPerfcondBottomPolShell(const Scalar_t lower, const Scalar_t upper,
      const Position position, const int l);

   /**
    * @brief Destructor
    */
   ~InsulatingTopPerfcondBottomPolShell() = default;

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

} // namespace Boundary
} // namespace LinearMap
} // namespace Chebyshev
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_INSULATINGTOPPERFCONDBOTTOMPOLSHELL_HPP
