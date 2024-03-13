/**
 * @file InsulatingSphereDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland insulating
 * sphere boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_INSULATINGSPHEREDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_INSULATINGSPHEREDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphereDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

/**
 * @brief Implementation of the full sphere Worland insulating sphere boundary
 * condition stencil
 */
class InsulatingSphereDiags
    : public QuICC::SparseSM::Worland::Stencil::InsulatingSphereDiags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   jacobi alpha
    * @param l       harmonic degree l
    */
   InsulatingSphereDiags(const Scalar_t alpha, const int l);

   /**
    * @brief Destructor
    */
   virtual ~InsulatingSphereDiags() = default;

   /**
    * @brief 1. subdiagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d_1(const ACoeff_t& n) const final;

   /**
    * @brief diagonal
    *
    * @param n Array of n indexes
    */
   ACoeff_t d0(const ACoeff_t& n) const final;

protected:
private:
};

} // namespace SphEnergy
} // namespace Stencil
} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_INSULATINGSPHEREDIAGS_HPP
