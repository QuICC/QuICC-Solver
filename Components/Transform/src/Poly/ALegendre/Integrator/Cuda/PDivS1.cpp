/**
 * @file PDivS1.cpp
 * @brief Source of the implementation of the associated Legendre P parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivS1.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/sin_1Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/InnerProduct.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<>
   PDivS1<CudaIALegendreOperatorTypes>::PDivS1()
   {
   }

   template<>
   PDivS1<CudaIALegendreOperatorTypes>::~PDivS1()
   {
   }

   template <>
   void PDivS1<CudaIALegendreOperatorTypes>::makeOperator(OpMatrix &op,
      const OpArray &igrid, const OpArray &iweights, const int i) const {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::sin_1Plm dplm;
      dplm.compute<MHDFloat>(op, nPoly, m, igrid, iweights, ev::Set());
   }

   template <>
   void PDivS1<CudaIALegendreOperatorTypes>::applyOperator(
      OpMatrixR rOut, const int i, const OpMatrixCR &in) const {
#if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
#endif // defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

#ifdef KOKKOS_ENABLE_CUDA
   template <>
   void PDivS1<CudaIALegendreOperatorTypes>::applyUnitOperator(
      const OpMatrixLZ &rOutView, const OpMatrixLZ &inView,
      const OpVectorI &scan, const int total) const {
      applyBlockOperator(
         this->mspSetup, this->vmOps, rOutView, inView, scan, total);
   }
#endif

template class PDivS1<CudaIALegendreOperatorTypes>;
} // namespace Integrator
}
}
}
}
