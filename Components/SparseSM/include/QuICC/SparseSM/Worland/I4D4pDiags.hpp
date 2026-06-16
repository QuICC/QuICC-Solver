/**
 * @file I4D4pDiags.hpp
 * @brief Interface to I4D4p diagonals for full sphere Worland I4D4p sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4D4PDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I4D4PDIAGS_HPP

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
 * @brief Implementation of the full sphere Worland I4D4p sparse operator.
 *
 * Band: subdiagonals d_4..d_1, main d0 (l+4 coupling).
 */
class I4D4pDiags : public IDiags
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
   I4D4pDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
      const int q);

   /**
    * @brief Destructor
    */
   virtual ~I4D4pDiags() = default;

   /**
    * @brief 4. subdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d_4(const ACoeff_t& n) const = 0;

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

protected:
private:
};

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_I4D4PDIAGS_HPP
