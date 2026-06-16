/**
 * @file I4laplD2pDiags.hpp
 * @brief Interface to I4laplD2p diagonals for full sphere Worland I4laplD2p
 * sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4LAPLD2PDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I4LAPLD2PDIAGS_HPP

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
 * @brief Implementation of the full sphere Worland I4laplD2p sparse operator.
 *
 * Band: subdiagonals d_3..d_1, main d0, superdiagonal d1 (l+2 coupling).
 */
class I4laplD2pDiags : public IDiags
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
   I4laplD2pDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
      const int q);

   /**
    * @brief Destructor
    */
   virtual ~I4laplD2pDiags() = default;

   /**
    * @brief 3. subdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d_3(const ACoeff_t& n) const = 0;

   /**
    * @brief 2. subdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d_2(const ACoeff_t& n) const = 0;

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

protected:
private:
};

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_I4LAPLD2PDIAGS_HPP
