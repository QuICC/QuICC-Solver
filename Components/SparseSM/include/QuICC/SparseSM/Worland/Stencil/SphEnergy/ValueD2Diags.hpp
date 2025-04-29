/**
 * @file ValueD2Diags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland ValueD2 boundary
 * condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_VALUED2DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_VALUED2DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/ValueD2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace SphEnergy {

/**
 * @brief Implementation of the full sphere Worland ValueD2 boundary condition
 * stencil
 */
class ValueD2Diags : public QuICC::SparseSM::Worland::Stencil::ValueD2Diags
{
public:
   /**
    * @brief Constructor
    *
    * @param alpha   jacobi alpha
    * @param l       harmonic degree l
    */
   ValueD2Diags(const Scalar_t alpha, const int l);

   /**
    * @brief Destructor
    */
   virtual ~ValueD2Diags() = default;

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

#endif // QUICC_SPARSESM_WORLAND_STENCIL_SPHENERGY_VALUED2DIAGS_HPP
