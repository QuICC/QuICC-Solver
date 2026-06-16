/**
 * @file I3laplDmDiags.hpp
 * @brief Interface to I3laplDm diagonals for full sphere Worland I3laplDm
 * sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I3LAPLDMDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I3LAPLDMDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/IDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

/**
 * @brief Implementation of the full sphere Worland I3laplDm sparse operator.
 *
 * Band: subdiagonal d_1, main d0, superdiagonals d1..d2 (l-1 coupling).
 */
class I3laplDmDiags : public IDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param dBeta   Jacobi beta = l + dBeta
    * @param l       Harmonic degree l
    * @param q       Truncation q (only consider rows - q equations)
    */
   I3laplDmDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
      const int q);

   /**
    * @brief Destructor
    */
   virtual ~I3laplDmDiags() = default;

   /**
    * @brief 1. subdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d_1(const ACoeff_t& n) const = 0;

   /**
    * @brief Main diagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d0(const ACoeff_t& n) const = 0;

   /**
    * @brief 1. superdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d1(const ACoeff_t& n) const = 0;

   /**
    * @brief 2. superdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d2(const ACoeff_t& n) const = 0;

protected:
private:
};

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_I3LAPLDMDIAGS_HPP
