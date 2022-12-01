/**
 * @file IALegendreIntegrator.cpp
 * @brief Source of the interface to a associated Legendre based integrator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   template<typename OpTypes>
   IALegendreIntegrator<OpTypes>::IALegendreIntegrator()
   {
   }

   template<typename OpTypes>
   IALegendreIntegrator<OpTypes>::~IALegendreIntegrator()
   {
   }


   template<typename OpTypes>
   void IALegendreIntegrator<OpTypes>::initOperators(const OpArray& igrid, const OpArray& iweights) const
   {
      // Initit specialized data for operators
      this->initSpecial();

      #if defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
      // reserve storage
      this->mOps.reserve(this->mspSetup->slowSize());

      // Loop over harmonic orders
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         // Build operator
         this->mOps.push_back(Matrix(igrid.size(), this->mspSetup->fastSize(i)));
         this->makeOperator(this->mOps.back(), igrid, iweights, i);
      }
      #elif defined QUICC_ALEGENDRE_INTGIMPL_OTF

         // Store grid and weights
         this->mGrid = igrid;
         this->mWeights = iweights;

      #endif //defined QUICC_ALEGENDRE_INTGIMPL_MATRIX
   }

   template<typename OpTypes>
   void IALegendreIntegrator<OpTypes>::applyOperators(OpMatrixZ& rOut, const OpMatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix("IALegendreIntegrator::applyOperators");

      // assert right sizes for input matrix
      assert(in.rows() == this->mspSetup->fwdSize());
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.cols() == this->mspSetup->blockSize());

      int start = 0;
      int inRows = this->mspSetup->fwdSize();
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int m = this->mspSetup->slow(i);
         int outRows = this->mspSetup->fast(this->mspSetup->fastSize(i)-1,i) - m + 1;
         assert(outRows == this->mspSetup->fastSize(i));

         this->applyOperator(rOut.block(0, start, outRows, cols), i, in.block(0, start, inRows, cols));

         start += cols;
      }
   }

   template<typename OpTypes>
   void IALegendreIntegrator<OpTypes>::initSpecial() const
   {
   }

   template<typename OpTypes>
   int IALegendreIntegrator<OpTypes>::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   template<typename OpTypes>
   int IALegendreIntegrator<OpTypes>::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   template<typename OpTypes>
   MHDFloat IALegendreIntegrator<OpTypes>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += IALegendreOperator::requiredStorage();

      // Storage for the operators
      for(auto it = this->mOps.cbegin(); it != this->mOps.cend(); ++it)
      {
         mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(it->size());
      }

      // Storage for grid and weights
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mGrid.size());
      mem += static_cast<MHDFloat>(Debug::MemorySize<MHDFloat>::BYTES)*(this->mWeights.size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   template class IALegendreIntegrator<IALegendreOperatorTypes>;
   template class IALegendreIntegrator<PIALegendreOperatorTypes>;
   template class IALegendreIntegrator<CudaIALegendreOperatorTypes>;
}
}
}
}
}
