/**
 * @file IWorlandProjector.cpp
 * @brief Source of the interface to a associated Worland based projector
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Projector/Kokkos/KokkosIWorlandProjector.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {


void KokkosIWorlandProjector::initOperators(const Internal::Array& igrid,
   const Internal::Array& iweights) const
{
   // Calc total number of columns for the big matrix
   auto total = 0;
   auto slowSize = this->mspSetup->slowSize();
   std::vector<int> scan(slowSize + 1, 0);

   for (int i = 0; i < slowSize; i++)
   {
      scan[i + 1] = scan[i] + this->mspSetup->fastSize(i);
   }
   total = scan[slowSize];
   // reserve storage on the device for horizontal layout left
   // It will not affect the cuda limit on the horizontal block dim.
   this->vmOps = OpMatrixL("vmops", igrid.size(), total);
   // Create host view
   auto vmOpsHost = Kokkos::create_mirror_view(this->vmOps);

   // Loop over harmonic orders
   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      // Build operator
      Matrix op;
      this->makeOperator(op, igrid, iweights, i);
      add_contribution_to_view_right(vmOpsHost, scan[i], op);
   }

   Kokkos::deep_copy(this->vmOps, vmOpsHost);
}

void KokkosIWorlandProjector::applyOperators(MatrixZ& rOut,
   const MatrixZ& in) const
{
   Profiler::RegionFixture<3> fix("PIWorlandProjector::applyOperators");

   // assert right sizes for input  matrix
   assert(in.cols() == this->mspSetup->blockSize());
   // assert right sizes for output matrix
   assert(rOut.rows() == this->outRows());
   assert(rOut.cols() == this->outCols());

   auto slowSize = this->mspSetup->slowSize();
   OpVectorI scan("outRows Scan", slowSize + 1);
   auto hostScan = Kokkos::create_mirror_view(scan);

   Profiler::RegionStart<4>("PIWorlandProjector::hostScan");

   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int cols = this->mspSetup->mult(i);
      hostScan[i + 1] = hostScan[i] + cols;
   }

   Profiler::RegionStop<4>("PIWorlandProjector::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);

   OpMatrixLZL inView("inView", in.rows(), in.cols());
   DeepCopyEigen(inView, in);

   OpMatrixLZ rOutView("rOutView", rOut.rows(), rOut.cols());

   Profiler::RegionStart<4>("PIWorlandProjector::applyUnitOperator");

   this->applyUnitOperator(rOutView, inView, scan, total);

   Profiler::RegionStop<4>("PIWorlandProjector::applyUnitOperator");

   DeepCopyEigen(rOut, rOutView);
}

void KokkosIWorlandProjector::defaultApplyOperator(OpMatrixR rOut, const int i,
   const OpMatrixCR& in) const
{
   rOut = this->mOps.at(i).transpose() * in;
}

int KokkosIWorlandProjector::outRows() const
{
   return this->mspSetup->fwdSize();
}

int KokkosIWorlandProjector::outCols() const
{
   return this->mspSetup->blockSize();
}

void KokkosIWorlandProjector::applyOperator(OpMatrixR rOut, const int i,
   const OpMatrixCR& in) const
{
   throw std::logic_error("P AL operator should never call this");
}

} // namespace Projector
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
