/**
 * @file IALegendreProjector.cpp
 * @brief Source of the interface to a associated Legendre based projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/Kokkos/KokkosIALegendreProjector.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Profiler/Interface.hpp"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   void KokkosIALegendreProjector::initSpecial() const {}

   int KokkosIALegendreProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int KokkosIALegendreProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   void KokkosIALegendreProjector::initOperators(const OpArray& igrid, const OpArray& iweights) const
   {
      // Initit specialized data for operators
      this->initSpecial();

      //Calc total number of columns for the big matrix
      auto total = 0;
      auto slowSize = this->mspSetup->slowSize();
      std::vector<int> scan(slowSize + 1, 0);

      for(int i = 0; i < slowSize; i++)
      {
         scan[i + 1] = scan[i] + this->mspSetup->fastSize(i);
      }
      total = scan[slowSize];
      // reserve storage on the device for horizontal layout left
      // It will not affect the cuda limit on the horizontal block dim.
      this->vmOps = OpMatrixL("vmops", igrid.size(), total);
      //Create host view
      auto vmOpsHost= Kokkos::create_mirror_view(this->vmOps);

      // Loop over harmonic orders
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         Matrix op;
         this->makeOperator(op, igrid, iweights, i);
         add_contribution_to_view_right(vmOpsHost, scan[i], op);
      }

      Kokkos::deep_copy(vmOps, vmOpsHost);
   }

   void KokkosIALegendreProjector::applyOperators(
      OpMatrixZ &rOut, const OpMatrixZ &in) const
   {
      Profiler::RegionFixture<3> fix("KokkosIALegendreProjector::applyOperators");

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->mspSetup->fwdSize());
      assert(rOut.cols() == this->mspSetup->blockSize());

      auto slowSize = this->mspSetup->slowSize();
      OpVectorI scan("outRows Scan", slowSize + 1);
      auto hostScan = Kokkos::create_mirror_view(scan);

      Profiler::RegionStart<4>("KokkosIALegendreProjector::hostScan");

      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int m = this->mspSetup->slow(i);
         int inRows =
            this->mspSetup->fast(this->mspSetup->fastSize(i) - 1, i) - m + 1;
         hostScan[i + 1] = hostScan[i] + inRows;
      }

      Profiler::RegionStop<4>("KokkosIALegendreProjector::hostScan");

      Kokkos::deep_copy(scan, hostScan);

      auto total = hostScan(slowSize);
      auto col_size = this->mspSetup->mult(0);

      // TODO: This needs to be done once in the beginning not need to be done
      // multiple times here.
      OpMatrixLZ inView("inView", total, col_size);
      DeepCopyEigen(inView, in, hostScan, col_size);

      OpMatrixLZ rOutView("rOutView", rOut.rows() * slowSize, col_size);

      Profiler::RegionStart<4>("KokkosIALegendreProjector::applyUnitOperator");

      this->applyUnitOperator(rOutView, inView, scan, total);

      Profiler::RegionStop<4>("KokkosIALegendreProjector::applyUnitOperator");

      DeepCopyEigen(rOut, rOutView, col_size);
   }

   void KokkosIALegendreProjector::applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const
   {
      throw std::logic_error("P AL operator should never call this");
   }

}
}
}
}
}
