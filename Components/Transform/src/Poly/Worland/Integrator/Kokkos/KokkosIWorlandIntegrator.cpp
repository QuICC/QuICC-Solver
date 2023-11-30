/**
 * @file IWorlandIntegrator.cpp
 * @brief Source of the interface to a Worland based integrator
 */

// System includes
//

// Project includes
#include "Profiler/Interface.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/Kokkos/KokkosIWorlandIntegrator.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

void KokkosIWorlandIntegrator::initOperators(const Internal::Array& igrid,
   const Internal::Array& iweights) const
{
   // Reserve storage for the operators
   auto total = 0;
   auto slowSize = this->mspSetup->slowSize();
   std::vector<Integer> scan(slowSize + 1, 0);

   for (int i = 0; i < slowSize; i++)
   {
      scan[i + 1] = scan[i] + this->mspSetup->fastSize(i);
   }
   total = scan[slowSize];

   // reserve storage on the device
   this->vmOps = OpMatrixL("vmops", total, igrid.size());
   // Create host view
   auto vmOpsHost = Kokkos::create_mirror_view(this->vmOps);

   // Loop over harmonic degrees
   for (int i = 0; i < slowSize; i++)
   {
      // Build operator
      Matrix op;
      this->makeOperator(op, igrid, iweights, i);
      auto transpose = op.transpose();
      add_contribution_to_view_left(vmOpsHost, scan[i], transpose);
   }

   Kokkos::deep_copy(vmOps, vmOpsHost);
}

void KokkosIWorlandIntegrator::applyOperators(MatrixZ& rOut,
   const MatrixZ& in) const
{
   Profiler::RegionFixture<3> fix("KokkosIWorlandIntegrator::applyOperators");

   // assert right sizes for input matrix
   assert(in.rows() == this->mspSetup->fwdSize());
   assert(in.cols() == this->mspSetup->blockSize());
   // assert right sizes for output matrix
   assert(rOut.cols() == this->mspSetup->blockSize());

   auto slowSize = this->mspSetup->slowSize();
   OpVectorI scan("outRows Scan", slowSize + 1);
   auto hostScan = Kokkos::create_mirror_view(scan);

   Profiler::RegionStart<4>("KokkosIWorlandIntegrator::hostScan");

   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int cols = this->mspSetup->mult(i);
      hostScan[i + 1] = hostScan[i] + cols;
   }

   Profiler::RegionStop<4>("KokkosIWorlandIntegrator::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);

   OpMatrixLZL inView("inView", in.rows(), in.cols());
   DeepCopyEigen(inView, in);

   OpMatrixLZ rOutView("rOutView", rOut.rows(), rOut.cols());

   Profiler::RegionStart<4>("KokkosIWorlandIntegrator::applyUnitOperator");

   this->applyUnitOperator(rOutView, inView, scan, total);

   Profiler::RegionStop<4>("KokkosIWorlandIntegrator::applyUnitOperator");

   DeepCopyEigen(rOut, rOutView);
}

int KokkosIWorlandIntegrator::outRows() const
{
   return this->mspSetup->fastSize(0);
}

int KokkosIWorlandIntegrator::outCols() const
{
   return this->mspSetup->blockSize();
}

void KokkosIWorlandIntegrator::applyOperator(OpMatrixR rOut, const int i,
   const OpMatrixCR& in) const
{
   throw std::logic_error("P Worland operator should never call this");
}

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
