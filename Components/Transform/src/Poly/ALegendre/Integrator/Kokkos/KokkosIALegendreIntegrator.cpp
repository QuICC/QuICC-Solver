/**
 * @file IALegendreIntegrator.cpp
 * @brief Source of the interface to a associated Legendre based integrator
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/Kokkos/KokkosIALegendreIntegrator.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

void KokkosIALegendreIntegrator::initSpecial() const {}

int KokkosIALegendreIntegrator::outRows() const
{
   return this->mspSetup->fastSize(0);
}

int KokkosIALegendreIntegrator::outCols() const
{
   return this->mspSetup->blockSize();
}

void KokkosIALegendreIntegrator::initOperators(const Internal::Array& igrid,
   const Internal::Array& iweights) const
{
   // Initit specialized data for operators
   this->initSpecial();

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

   // Loop over harmonic orders
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

void KokkosIALegendreIntegrator::applyOperators(MatrixZ& rOut,
   const MatrixZ& in) const
{
   Profiler::RegionFixture<3> fix("KokkosIALegendreIntegrator::applyOperators");

   // assert right sizes for input matrix
   assert(in.rows() == this->mspSetup->fwdSize());
   assert(in.cols() == this->mspSetup->blockSize());
   // assert right sizes for output matrix
   assert(rOut.cols() == this->mspSetup->blockSize());

   auto slowSize = this->mspSetup->slowSize();
   OpVectorI scan("outRows Scan", slowSize + 1);
   auto hostScan = Kokkos::create_mirror_view(scan);

   Profiler::RegionStart<4>("KokkosIALegendreIntegrator::hostScan");

   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int m = this->mspSetup->slow(i);
      int outRows =
         this->mspSetup->fast(this->mspSetup->fastSize(i) - 1, i) - m + 1;
      hostScan[i + 1] = hostScan[i] + outRows;
   }

   Profiler::RegionStop<4>("KokkosIALegendreIntegrator::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);

   // TODO: The following deep copy of the in data needs to be done once.
   // It is called multiple times here unnecessarily
   OpMatrixLZ inView("inView", in.rows(), in.cols());
   DeepCopyEigen(inView, in);

   // TODO: Vertical "in" matrix option in case we decide to go this way.
   /* OpMatrixLZ inView("inView", in.rows() * slowSize, col_size);
   DeepCopyEigen(inView, in, col_size); */

   auto col_size = this->mspSetup->mult(0);
   OpMatrixLZ rOutView("rOutView", total, col_size);

   Profiler::RegionStart<4>("KokkosIALegendreIntegrator::applyUnitOperator");

   this->applyUnitOperator(rOutView, inView, scan, total);

   Profiler::RegionStop<4>("KokkosIALegendreIntegrator::applyUnitOperator");

   DeepCopyEigen(rOut, rOutView, hostScan, col_size);
}

void KokkosIALegendreIntegrator::applyOperator(OpMatrixR rOut, const int i,
   const OpMatrixCR& in) const
{
   throw std::logic_error("P AL operator should never call this");
}

} // namespace Integrator
} // namespace ALegendre
} // namespace Poly
} // namespace Transform
} // namespace QuICC
