/**
 * @file D1DivS1Dp.cpp
 * @brief Source of the implementation of the associated Legendre 1/S D_phi
 * projector
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/D1DivS1Dp.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Polynomial/ALegendre/dsin_1Plm.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

void D1DivS1Dp<kokkos_t>::applyUnitOperator(const OpMatrixLZ& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
   constantMultiplyMatrix<1>(this->mspSetup, scan, inView);
   D1DivS1<kokkos_t>::applyUnitOperator(rOutView, inView, scan, total);
}

} // namespace Projector
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
