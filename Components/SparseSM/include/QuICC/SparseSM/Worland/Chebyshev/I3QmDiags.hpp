/**
 * @file I3QmDiags.hpp
 * @brief Interface to I3Qm diagonals for full sphere Worland I3Qm sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3QMDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3QMDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/I3QmDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

/**
 * @brief Implementation of the full sphere Worland I3Qm sparse operator
 */
class I3QmDiags : public QuICC::SparseSM::Worland::I3QmDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param l       Harmonic degree
    * @param q       Truncation q
    */
   I3QmDiags(const Scalar_t alpha, const int l, const int q);

   /**
    * @brief Destructor
    */
   virtual ~I3QmDiags() = default;

   /**
    * @brief 2. subdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d_2(const ACoeff_t& n) const final;

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

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3QMDIAGS_HPP
