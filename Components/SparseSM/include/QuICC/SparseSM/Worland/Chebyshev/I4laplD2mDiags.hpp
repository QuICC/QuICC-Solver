/**
 * @file I4laplD2mDiags.hpp
 * @brief Interface to I4laplD2m diagonals for full sphere Worland I4laplD2m
 * sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4LAPLD2MDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4LAPLD2MDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/I4laplD2mDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

/**
 * @brief Implementation of the full sphere Worland I4laplD2m sparse operator
 *
 * Velocity-only spheroid LHS operator: no l=0 patch (l>=1, Id at l=0 in the
 * backend), no correctQ1/zeroLast (TRUNCATE_QI=OFF; q>0 unsupported).
 */
class I4laplD2mDiags : public QuICC::SparseSM::Worland::I4laplD2mDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param l       Harmonic degree
    * @param q       Truncation q
    */
   I4laplD2mDiags(const Scalar_t alpha, const int l, const int q);

   /**
    * @brief Destructor
    */
   virtual ~I4laplD2mDiags() = default;

   /**
    * @brief 1. subdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d_1(const ACoeff_t& n) const final;

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

protected:
private:
};

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I4LAPLD2MDIAGS_HPP
