/**
 * @file PP.cpp
 * @brief Source of the implementation of the associated Legendre PP projector
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Projector/PP.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   template <> PP<PIALegendreOperatorTypes>::PP() {}

   template <> PP<PIALegendreOperatorTypes>::~PP() {}

   template<>
   void PP<PIALegendreOperatorTypes>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, OpArray(), ev::Set());
   }


   template<>
   void PP<PIALegendreOperatorTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

   template <typename OpTypes>
   void PP<OpTypes>::applyOperator(
      OpMatrixR rOut, const int i, const OpMatrixCR &in) const {
#if defined QUICC_ALEGENDRE_PROJIMPL_MATRIX
      OpMatrixL mOpsView(
         "mOpsView", this->mOps.at(i).rows(), this->mOps.at(i).cols());
      DeepCopyEigen(mOpsView, this->mOps.at(i));

      OpMatrixLZ inView("inView", in.rows(), in.cols());
      DeepCopyEigen(inView, in);

      OpMatrixLZ rOutView("rOutView", rOut.rows(), rOut.cols());
      denseMatrixMultiply("T", "N", mOpsView, inView, rOutView);

      DeepCopyEigen(rOut, rOutView);
#endif
   }


   //Requires rOut to be a vertical matrix instead of horizontal.
   //Better coalescing achieved with very good occupancy.
   //However requires specialized copy of rout into its original format
   template <typename R, typename T, typename V> struct ApplyUnitOperator {

      ApplyUnitOperator(R rv, R iv, T vo, V s, int tl, int ir)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), slowSize(tl),
            outRows(ir) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = row / outRows;
         auto local_row = row % outRows;

         auto in_rows_start = scan(index);
         auto in_rows_size = scan(index + 1) - in_rows_start;

         for(int j = 0; j < in_rows_size; j++)
         {
            auto inner_index = j + in_rows_start;
            rOutView(row, col) +=
               Ops(local_row, inner_index) * inView(inner_index, col);
         }
      }

      R rOutView;
      R inView;
      T Ops;
      V scan;
      int slowSize;
      int outRows;
   };


   template <typename OpTypes>
   void PP<OpTypes>::applyUnitOperator(const OpMatrixLZ &rOutView,
      const OpMatrixLZ &inView, const OpVectorI &scan, const int total) const {
      auto outRows = this->mspSetup->fwdSize();
      auto slowSize = this->mspSetup->slowSize();

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
            {0, 0}, {rOutView.extent(0), rOutView.extent(1)}),
         ApplyUnitOperator(
            rOutView, inView, this->vmOps, scan, slowSize, outRows));
   }

   template class PP<PIALegendreOperatorTypes>;
}
}
}
}
}
