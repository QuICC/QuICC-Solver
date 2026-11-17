/**
 * @file ComplexFourierTransform.cpp
 * @brief Source of the implementation of the FFTW transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/ComplexFourierTransform.hpp"

#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D2.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D3.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/D1.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Df1Lapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Ds1Lapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Lapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/Mean.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/P.hpp"

#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P_Clean.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D1_P.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D2.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/D1.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Df1InvLapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/InvLapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Lapl2D.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/Mean.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/P0.hpp"
#include "QuICC/Transform/Forward/Pm.hpp"
#include "QuICC/Transform/Forward/D1.hpp"
#include "QuICC/Transform/Forward/D1ZP0.hpp"
#include "QuICC/Transform/Forward/D2.hpp"
#include "QuICC/Transform/Forward/Laplh.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/DfOverlaplh.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/P0.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/D3.hpp"
#include "QuICC/Transform/Backward/DfLaplh.hpp"
#include "QuICC/Transform/Backward/DsLaplh.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"

namespace QuICC {

namespace Transform {

   void ComplexFourierTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void ComplexFourierTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array ComplexFourierTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void ComplexFourierTransform::init(ComplexFourierTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise FFTW plans
      this->initOperators();
   }

   void ComplexFourierTransform::initOperators()
   {
      using namespace ::QuICC::Transform::Fft::Fourier;
      #ifdef QUICC_HAS_CUDA_BACKEND
         #ifdef QUICC_USE_VKFFT
            using backend_t = viewGpuVkFFT_t;
         #else
            #ifdef QUICC_USE_CUFFT
               using backend_t = viewGpu_t;
            #endif
         #endif
      #else
         using backend_t = base_t;
      #endif
      
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::P<backend_t>>(Backward::P::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Mean<backend_t>>(Backward::P0::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D1<backend_t>>(Backward::D1::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D2<backend_t>>(Backward::D2::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D3<backend_t>>(Backward::D3::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Df1Lapl2D<backend_t>>(Backward::DfLaplh::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Ds1Lapl2D<backend_t>>(Backward::DsLaplh::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Lapl2D<backend_t>>(Backward::Laplh::id());

         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::P<backend_t>>(Forward::P::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::P_Clean<backend_t>>(Forward::Pm::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D1_P<backend_t>>(Forward::D1ZP0::id());

         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D1<backend_t>>(Forward::D1::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D2<backend_t>>(Forward::D2::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Df1InvLapl2D<backend_t>>(Forward::DfOverlaplh::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::InvLapl2D<backend_t>>(Forward::Overlaplh::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Lapl2D<backend_t>>(Forward::Laplh::id());
         this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Mean<backend_t>>(Forward::P0::id());
   }

   void ComplexFourierTransform::forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void ComplexFourierTransform::backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   MHDFloat ComplexFourierTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void ComplexFourierTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRACOMPLEXFFT, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
