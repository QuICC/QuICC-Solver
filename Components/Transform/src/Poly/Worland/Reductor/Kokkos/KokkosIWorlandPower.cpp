/**
 * @file IWorlandPower.cpp
 * @brief Source of the interface to a Worland based reduction operator (e.g.
 * energy)
 */

// System includes
//

// project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandPower.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

KokkosIWorlandPower::KokkosIWorlandPower(const int shift) : mcShift(shift) {}

KokkosIWorlandPower::~KokkosIWorlandPower() {}

void KokkosIWorlandPower::initOperators(const Internal::Array& icompgrid,
   const Internal::Array& icompweights) const
{
   // Energy calculation requires a different quadrature
   Internal::Array igrid, iweights;
   this->computePowerQuadrature(igrid, iweights, icompgrid.size());

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
   this->vmOps = OpMatrixL("vmOps", igrid.size(), total);
   this->vmEOps = OpMatrixL("vmEOps", total, igrid.size());
   // Create host view
   auto vmOpsHost = Kokkos::create_mirror_view(this->vmOps);
   auto vmEOpsHost = Kokkos::create_mirror_view(this->vmEOps);

   // Loop over harmonic orders
   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      // Build operator
      Matrix op, eop;
      this->makeOperator(op, eop, igrid, iweights, i);
      auto transpose = eop.transpose();
      add_contribution_to_view_right(vmOpsHost, scan[i], op);
      add_contribution_to_view_left(vmEOpsHost, scan[i], transpose);
   }

   Kokkos::deep_copy(this->vmOps, vmOpsHost);
   Kokkos::deep_copy(this->vmEOps, vmEOpsHost);
}

void KokkosIWorlandPower::computePowerQuadrature(Internal::Array& igrid,
   Internal::Array& iweights, const int gSize) const
{
   int nrgSize = gSize + 2 * this->mcShift;

   Polynomial::Quadrature::WorlandSphEnergyRule wquad;
   wquad.computeQuadrature(igrid, iweights, nrgSize);
}

void KokkosIWorlandPower::applyOperators(Matrix& rOut, const MatrixZ& in) const
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

   for (int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      int cols = this->mspSetup->mult(i);
      hostScan[i + 1] = hostScan[i] + cols;
   }

   Profiler::RegionStop<4>("KokkosIWorlandProjector::hostScan");

   Kokkos::deep_copy(scan, hostScan);

   auto total = hostScan(slowSize);
   OpMatrixLZL inView("inView", in.rows(), in.cols());
   DeepCopyEigen(inView, in);
   OpMatrixL rOutView("rOutView", rOut.rows(), rOut.cols());

   Profiler::RegionStart<4>("KokkosIWorlandProjector::applyUnitOperator");
   this->applyUnitOperator(rOutView, inView, scan, total);
   Profiler::RegionStop<4>("KokkosIWorlandProjector::applyUnitOperator");

   DeepCopyEigen(rOut, rOutView);
}

void KokkosIWorlandPower::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
{
   throw std::logic_error("Unused interface");
}

void KokkosIWorlandPower::defaultApplyUnitOperator(const OpMatrixL& rOutView,
   const OpMatrixLZ& inView, const OpVectorI& scan, const int total) const
{
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
   OpMatrixLZ temp("Temp View", this->vmOps.extent(0), inView.extent(1));
   applyBlockOperator<4>(this->mspSetup, this->vmOps, temp, inView, scan,
      total);
   applyBlockOperator<3>(this->mspSetup, this->vmEOps, rOutView, temp, scan,
      total, Abs2Complex());
#endif
}

int KokkosIWorlandPower::outRows() const
{
   return this->mspSetup->fastSize(0);
}

int KokkosIWorlandPower::outCols() const
{
   return this->mspSetup->blockSize();
}

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC
