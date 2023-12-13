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
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/P.hpp"

// Project includes
//
#include "QuICC/Polynomial/ALegendre/Plm.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/Set.hpp"
#include "QuICC/Polynomial/ALegendre/Evaluator/OuterProduct.hpp"

#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void P<kokkos_t>::makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      int m = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1 ;

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::ALegendre::Evaluator;
      Polynomial::ALegendre::Plm plm;
      plm.compute<MHDFloat>(op, nPoly, m, igrid, Internal::Array(), ev::Set());
   }


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


   void P<kokkos_t>::applyUnitOperator(const OpMatrixLZ &rOutView,
      const OpMatrixLZ &inView, const OpVectorI &scan, const int total) const {
      auto outRows = this->mspSetup->fwdSize();
      auto slowSize = this->mspSetup->slowSize();

      Kokkos::parallel_for("Apply Kokkos Operator",
         Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
            {0, 0}, {rOutView.extent(0), rOutView.extent(1)}),
         ApplyUnitOperator(
            rOutView, inView, this->vmOps, scan, slowSize, outRows));
   }

}
}
}
}
}
