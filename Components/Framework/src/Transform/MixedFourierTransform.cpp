/**
 * @file MixedFourierTransform.cpp
 * @brief Source of the implementation of the FFTW mixed transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/MixedFourierTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D1.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D2.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D3.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/D4.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/P.hpp"

#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1_P.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D1.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/D2.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/D1.hpp"
#include "QuICC/Transform/Forward/D1ZP0.hpp"
#include "QuICC/Transform/Forward/D2.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/D3.hpp"

#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   MixedFourierTransform::MixedFourierTransform()
   {
   }

   MixedFourierTransform::~MixedFourierTransform()
   {
   }

   void MixedFourierTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void MixedFourierTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array MixedFourierTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void MixedFourierTransform::init(MixedFourierTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise FFTW plans
      this->initOperators();
   }

   void MixedFourierTransform::initOperators()
   {
      using namespace QuICC::Transform::Fft::Fourier;
      this->mImpl.addOperator<Fft::Fourier::Mixed::Projector::P<base_t>>(Backward::P::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Projector::D1<base_t>>(Backward::D1::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Projector::D2<base_t>>(Backward::D2::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Projector::D3<base_t>>(Backward::D3::id());

      this->mImpl.addOperator<Fft::Fourier::Mixed::Integrator::P<base_t>>(Forward::P::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Integrator::D1<base_t>>(Forward::D1::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Integrator::D1_P<base_t>>(Forward::D1ZP0::id());
      this->mImpl.addOperator<Fft::Fourier::Mixed::Integrator::D2<base_t>>(Forward::D2::id());
   }

   void MixedFourierTransform::forward(MatrixZ& rOut, const Matrix& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void MixedFourierTransform::backward(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("MixedFourierTransform::backward");
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void MixedFourierTransform::forward(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::backward(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::reduce(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void MixedFourierTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat MixedFourierTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void MixedFourierTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRAMIXEDFFT, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
