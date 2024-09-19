/**
 * @file I2LaplDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland I2Lapl sparse
 * operator
 */

#ifndef QUICC_SPARSESM_WORLAND_SPHENERGY_I2LAPLDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_SPHENERGY_I2LAPLDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/I2LaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

/**
 * @brief Implementation of the full sphere Worland I2Lapl sparse operator
 */
class I2LaplDiags : public QuICC::SparseSM::Worland::I2LaplDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   Jacobi alpha
    * @param l       Harmonic degree
    * @param q       Truncation q
    */
   I2LaplDiags(const Scalar_t alpha, const int l, const int q);

   /**
    * @brief Destructor
    */
   virtual ~I2LaplDiags() = default;

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

protected:
private:
};

} // namespace SphEnergy
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_SPHENERGY_I2LAPLDIAGS_HPP
