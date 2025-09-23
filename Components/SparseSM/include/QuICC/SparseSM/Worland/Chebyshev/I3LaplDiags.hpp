/**
 * @file I3LaplDiags.hpp
 * @brief Interface to I3Lapl diagonals for full sphere Worland I3Lapl sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3LAPLDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3LAPLDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I3Diags.hpp"
#include "QuICC/SparseSM/Worland/I3LaplDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

/**
 * @brief Implementation of the full sphere Worland I3Lapl sparse operator
 */
class I3LaplDiags : public QuICC::SparseSM::Worland::I3LaplDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param l       Harmonic degree
    * @param q       Truncation q
    */
   I3LaplDiags(const Scalar_t alpha, const int l, const int q);

   /**
    * @brief Destructor
    */
   virtual ~I3LaplDiags() = default;

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

   /**
    * @brief 4. superdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d4(const ACoeff_t& n) const final;

protected:
   /**
    * @brief Correct diagonal k for q = 2 truncation
    *
    * @param val  Diagonal values to correct
    * @param n    Array of n indexes
    * @param k    Diagonal k
    */
   void correctQ2(ACoeff_t& val, const ACoeff_t& n, const int k) const;

private:
   /**
    * @brief I3 diagonals for correction
    */
   I3Diags mI3;
};

} // namespace Chebyshev
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I3LAPLDIAGS_HPP
