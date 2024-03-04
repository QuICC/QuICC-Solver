/**
 * @file LlDivS1Dp.cpp
 * @brief Source of the parallel implementation of the associated Legendre
 * l(l+1)/Sin D_phi projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/LlDivS1Dp.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

void LlDivS1Dp<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid,
   const Internal::Array& iweights, const int i) const
{
   DivS1<kokkos_t>::makeOperator(op, igrid, iweights, i);
   op = op * (static_cast<MHDFloat>(this->mspSetup->slow(i)) *
                this->mLl.bottomRows(op.cols()))
                .asDiagonal();
}

void LlDivS1Dp<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);

   DataType constant;
   constant.real() = Math::cI.real();
   constant.imag() = Math::cI.imag();

   constantMultiplyMatrix(constant, rOutView);
}

void LlDivS1Dp<kokkos_t>::initSpecial() const
{
   // Initialise storage for l(l+1) factor
   this->mLl = Array::LinSpaced(this->mspSetup->specSize(), 0,
      this->mspSetup->specSize() - 1);
   this->mLl = this->mLl.array() * (this->mLl.array() + 1.0);
}

} // namespace Projector
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
