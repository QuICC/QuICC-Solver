/**
 * @file I4D4mDiags.hpp
 * @brief Interface to I4D4m diagonals for full sphere Worland I4D4m sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D4MDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D4MDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/I4D4mDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

/**
 * @brief Implementation of the full sphere Worland I4D4m sparse operator
 *
 * Velocity-only spheroid LHS operator: no l=0 patch (l>=1, Id at l=0 in the
 * backend), no correctQ1/zeroLast (TRUNCATE_QI=OFF; q>0 unsupported).
 */
class I4D4mDiags : public QuICC::SparseSM::Worland::I4D4mDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param l       Harmonic degree
    * @param q       Truncation q
    */
   I4D4mDiags(const Scalar_t alpha, const int l, const int q);

   /**
    * @brief Destructor
    */
   virtual ~I4D4mDiags() = default;

   /**
    * @brief Main diagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d0(const ACoeff_t& n) const final;

   /**
    * @brief 1. superdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d1(const ACoeff_t& n) const final;

   /**
    * @brief 2. superdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d2(const ACoeff_t& n) const final;

   /**
    * @brief 3. superdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d3(const ACoeff_t& n) const final;

   /**
    * @brief 4. superdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d4(const ACoeff_t& n) const final;

protected:
private:
};

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4D4MDIAGS_HPP
