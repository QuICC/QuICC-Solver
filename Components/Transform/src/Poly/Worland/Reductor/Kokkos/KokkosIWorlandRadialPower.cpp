/**
 * @file KokkosIWorlandRadialPower.cpp
 * @brief Source of the interface to a Worland radial grid power operator (e.g.
 * energy)
 */

// System includes
//

// Project includes
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandRadialPower.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

void KokkosIWorlandRadialPower::initOperators(const Internal::Array& igrid,
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

void KokkosIWorlandRadialPower::applyOperators(Matrix& rOut,
   const MatrixZ& in) const
{
   Profiler::RegionFixture<3> fix("KokkosIWorlandProjector::applyOperators");

   // assert right sizes for input  matrix
   assert(in.cols() == this->mspSetup->blockSize());
   // assert right sizes for output matrix
   assert(rOut.rows() == this->outRows());
   assert(rOut.cols() == this->outCols());

   auto slowSize = this->mspSetup->slowSize();
   OpVectorI scan("outRows Scan", slowSize + 1);
   auto hostScan = Kokkos::create_mirror_view(scan);

   Profiler::RegionStart<4>("KokkosIWorlandProjector::hostScan");

   // Changing col size or row size? If change then it behaves as a worland
   // projector
   int change = this->mspSetup->mult(slowSize - 1) - this->mspSetup->mult(0);
   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int val = 0;
      if (change)
         val = this->mspSetup->mult(i);
      else
         val = this->mspSetup->fastSize(i);

      hostScan[i + 1] = hostScan[i] + val;
   }

   Profiler::RegionStop<4>("KokkosIWorlandProjector::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);

   if (change)
   {
      OpMatrixLZL inView("inView", in.rows(), in.cols());
      DeepCopyEigen(inView, in);
      OpMatrixL rOutView("rOutView", rOut.rows(), rOut.cols());

      Profiler::RegionStart<4>("KokkosIWorlandProjector::applyUnitOperator");
      this->applyUnitOperator(rOutView, inView, scan, total);
      Profiler::RegionStop<4>("KokkosIWorlandProjector::applyUnitOperator");

      DeepCopyEigen(rOut, rOutView);
   }
   else
   {
      auto col_size = this->mspSetup->mult(0);

      OpMatrixLZ inView("inView", total, col_size);
      DeepCopyEigen(inView, in, hostScan, col_size);
      OpMatrixL rOutView("rOutView", rOut.rows() * slowSize, col_size);

      Profiler::RegionStart<4>("KokkosIWorlandProjector::applyUnitOperator");
      this->applyUnitOperator(rOutView, inView, scan, total);
      Profiler::RegionStop<4>("KokkosIWorlandProjector::applyUnitOperator");

      DeepCopyEigen(rOut, rOutView, col_size);
   }
}

void KokkosIWorlandRadialPower::applyOperators(MatrixZ& rOut,
   const MatrixZ& in) const
{
   throw std::logic_error("Unused interface");
}

int KokkosIWorlandRadialPower::outRows() const
{
   return this->mspSetup->fwdSize();
}

int KokkosIWorlandRadialPower::outCols() const
{
   return this->mspSetup->blockSize();
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
