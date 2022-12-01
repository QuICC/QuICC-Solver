/**
 * @file PP.cpp
 * @brief Source of the implementation of the associated Legendre P parallel integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PP.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/InnerProduct.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
/* #include <type_traits> */

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<typename OpTypes>
   PP<OpTypes>::PP()
   {
   }

   template<typename OpTypes>
   PP<OpTypes>::~PP()
   {
   }

   template<typename OpTypes>
   void PP<OpTypes>::makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, iweights, ev::Set());
   }

   template <typename OpTypes>
   void PP<OpTypes>::applyOperator(
      OpMatrixR rOut, const int i, const OpMatrixCR &in) const {
#if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
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

   template<>
   void PP<PIALegendreOperatorTypes>::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      rOut = this->mOps.at(i).transpose() * in;
      #endif //defined QUICC_ALEGENDRE_INTGIMPL_OTF
   }

   //Operator working on horizontal rout layout
   /* template<typename R, typename T, typename V>
   struct ApplyOperator {

      ApplyOperator(R rv, R iv, T vo, V s, int tl, int ir, int cs)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), total(tl), inRows(ir),
            col_size(cs) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int thid) const {
         auto local_row = thid % total;
         auto index = binary_search_range(scan, local_row);

         auto row = local_row - scan(index);
         auto col = thid / total + index * col_size;

         for(int j = 0; j < inRows; j++)
         { rOutView(row, col) += Ops(local_row, j) * inView(j, col); }
      }


      R rOutView;
      R inView;
      T Ops;
      V scan;
      int total;
      int inRows;
      int col_size;
   }; */

   //Requires rOut to be a vertical matrix instead of horizontal.
   //Better coalescing achieved with very good occupancy.
   //However requires specialized copy of rout into its original format
   template <typename R, typename T, typename V>
   struct ApplyUnitOperator{

      ApplyUnitOperator(R rv, R iv, T vo, V s, int ir, int cs)
          : rOutView(rv), inView(iv), Ops(vo), scan(s), inRows(ir),
            col_size(cs) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int row, const int col) const {
         auto index = binary_search_range(scan, row);
         auto colin = col + index * col_size;

         for(int j = 0; j < inRows; j++)
         { rOutView(row, col) += Ops(row, j) * inView(j, colin); }
      }

      R rOutView;
      R inView;
      T Ops;
      V scan;
      int inRows;
      int col_size;
   };


   template <typename OpTypes>
   void PP<OpTypes>::applyUnitOperator(const OpMatrixLZ &rOutView,
      const OpMatrixLZ &inView, const OpVectorI &scan, const int total) const {
      auto inRows = this->mspSetup->fwdSize();
      auto col_size = this->mspSetup->mult(0);

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
            {0, 0}, {rOutView.extent(0), rOutView.extent(1)}),
         ApplyUnitOperator(
            rOutView, inView, this->vmOps, scan, inRows, col_size));
   }

   template class PP<PIALegendreOperatorTypes>;
}
}
}
}
}
