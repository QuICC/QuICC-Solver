/**
 * @file PIALegendreIntegrator.cpp
 * @brief Source of the interface to a associated Legendre based integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/PIALegendreIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

#include "QuICC/Transform/Poly/KokkosUtils.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<typename OpTypes>
   PIALegendreIntegrator<OpTypes>::PIALegendreIntegrator()
   {
   }

   template<typename OpTypes>
   PIALegendreIntegrator<OpTypes>::~PIALegendreIntegrator()
   {
   }

   template <typename OpTypes>
   void PIALegendreIntegrator<OpTypes>::initSpecial() const {}

   template<typename OpTypes>
   int PIALegendreIntegrator<OpTypes>::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   template<typename OpTypes>
   int PIALegendreIntegrator<OpTypes>::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   template<typename OpTypes>
   void PIALegendreIntegrator<OpTypes>::initOperators(const OpArray& igrid, const OpArray& iweights) const
   {
      // Initit specialized data for operators
      this->initSpecial();

      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      auto total = 0;
      auto slowSize = this->mspSetup->slowSize();
      std::vector<Integer> scan(slowSize + 1, 0);

      for(int i = 0; i < slowSize; i++)
      {
          scan[i + 1] = scan[i] + this->mspSetup->fastSize(i);
      }
      total = scan[slowSize];
      // reserve storage on the device
      this->vmOps = OpMatrixL("vmops", total, igrid.size());
      //Create host view
      auto vmOpsHost= Kokkos::create_mirror_view(this->vmOps);

      // Loop over harmonic orders
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         Matrix op;
         this->makeOperator(op, igrid, iweights, i);
         auto transpose = op.transpose();
         add_contribution_to_view_left(vmOpsHost, scan[i], transpose);
      }

      Kokkos::deep_copy(vmOps, vmOpsHost);

      #elif defined QUICC_ALEGENDRE_INTGIMPL_OTF

         // Store grid and weights
         this->mGrid = igrid;
         this->mWeights = iweights;

      #endif //defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
   }

   template<typename OpTypes>
   void PIALegendreIntegrator<OpTypes>::applyOperators(OpMatrixZ& rOut, const OpMatrixZ& in) const
   {
      // assert right sizes for input matrix
      assert(in.rows() == this->mspSetup->fwdSize());
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.cols() == this->mspSetup->blockSize());

      auto slowSize = this->mspSetup->slowSize();
      OpVectorI scan("outRows Scan", slowSize + 1);
      auto hostScan = Kokkos::create_mirror_view(scan);

      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int m = this->mspSetup->slow(i);
         int outRows = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1;
         hostScan[i + 1] = hostScan[i] + outRows;
      }

      Kokkos::deep_copy(scan, hostScan);

#if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      auto total = hostScan(slowSize);

      //TODO: The following deep copy of the in data needs to be done once.
      //It is called multiple times here unnecessarily
      OpMatrixLZ inView("inView", in.rows(), in.cols());
      DeepCopyEigen(inView, in);

      //TODO: Vertical "in" matrix option in case we decide to go this way.
      /* OpMatrixLZ inView("inView", in.rows() * slowSize, col_size);
      DeepCopyEigen(inView, in, col_size); */

      auto col_size = this->mspSetup->mult(0);
      OpMatrixLZ rOutView("rOutView", total, col_size);

      this->applyUnitOperator(rOutView, inView, scan, total);

      DeepCopyEigen(rOut, rOutView, hostScan, col_size);
#endif
   }

   template class PIALegendreIntegrator<PIALegendreOperatorTypes>;
   template class PIALegendreIntegrator<CudaIALegendreOperatorTypes>;
}
}
}
}
}
