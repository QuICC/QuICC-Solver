/** 
 * @file ComplexFourierTransform.cpp
 * @brief Source of the implementation of the FFTW transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/ComplexFourierTransform.hpp"

// Project includes
//

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
#include "QuICC/Transform/Forward/Laplh_1.hpp"
#include "QuICC/Transform/Forward/DfLaplh_1.hpp"

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

   ComplexFourierTransform::ComplexFourierTransform()
   {
   }

   ComplexFourierTransform::~ComplexFourierTransform()
   {
   }

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
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Mean>(Backward::P0::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D2>(Backward::D2::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::D3>(Backward::D3::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Df1Lapl2D>(Backward::DfLaplh::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Ds1Lapl2D>(Backward::DsLaplh::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Projector::Lapl2D>(Backward::Laplh::id());

      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::P_Clean>(Forward::Pm::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D1_P>(Forward::D1ZP0::id());

      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D1>(Forward::D1::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::D2>(Forward::D2::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Df1InvLapl2D>(Forward::DfLaplh_1::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::InvLapl2D>(Forward::Laplh_1::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Lapl2D>(Forward::Laplh::id());
      this->mImpl.addOperator<Fft::Fourier::Complex::Integrator::Mean>(Forward::P0::id());
   }

   void ComplexFourierTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void ComplexFourierTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void ComplexFourierTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::reduce(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ComplexFourierTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
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
