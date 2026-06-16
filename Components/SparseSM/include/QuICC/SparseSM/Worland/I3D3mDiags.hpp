/**
 * @file I3D3mDiags.hpp
 * @brief Interface to I3D3m diagonals for full sphere Worland I3D3m sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I3D3MDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I3D3MDIAGS_HPP

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
 * @brief Implementation of the full sphere Worland I3D3m sparse operator.
 *
 * Band: main d0, superdiagonals d1..d3 (l-3 coupling).
 */
class I3D3mDiags : public IDiags
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
   I3D3mDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
      const int q);

   /**
    * @brief Destructor
    */
   virtual ~I3D3mDiags() = default;

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

   /**
    * @brief 3. superdiagonal
    *
    * @param n Array of n indexes
    */
   virtual ACoeff_t d3(const ACoeff_t& n) const = 0;

protected:
private:
};

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_I3D3MDIAGS_HPP
