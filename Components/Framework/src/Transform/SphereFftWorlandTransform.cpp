/**
 * file SphereFftWorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/SphereFftWorlandTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Fft/Worland/Projector/P.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/D1.hpp"
//#include "QuICC/Transform/Fft/Worland/Projector/D1R1.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/SphLapl.hpp"

#include "QuICC/Transform/Fft/Worland/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/P_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/R1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I4DivR1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I4DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2DivR1_Zero.hpp"
#include "QuICC/Transform/Fft/Worland/Integrator/I2DivR1D1R1_Zero.hpp"

#include "QuICC/Transform/Fft/Worland/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/Energy.hpp"

#include "QuICC/Transform/Fft/Worland/Reductor/PowerSLaplR2.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/Power.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Slapl.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I2P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/Q.hpp"
#include "QuICC/Transform/Forward/S.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/EnergySlaplR2.hpp"

#include "QuICC/Transform/Reductor/Power.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Reductor/PowerSlaplR2.hpp"

namespace QuICC {

namespace Transform {

   SphereFftWorlandTransform::SphereFftWorlandTransform()
   {
   }

   SphereFftWorlandTransform::~SphereFftWorlandTransform()
   {
   }

   void SphereFftWorlandTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void SphereFftWorlandTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array SphereFftWorlandTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void SphereFftWorlandTransform::init(SphereFftWorlandTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the ID to operator mapping
      this->initOperators();
   }

   void SphereFftWorlandTransform::initOperators()
   {
      using backend_t = Fft::Worland::base_t;

      // Create projectors
      this->mImpl.addOperator<Fft::Worland::Projector::P<backend_t>>(Backward::P::id());
      this->mImpl.addOperator<Fft::Worland::Projector::DivR1_Zero<backend_t>>(Backward::Overr1::id());
      this->mImpl.addOperator<Fft::Worland::Projector::D1<backend_t>>(Backward::D1::id());
//      this->mImpl.addOperator<Fft::Worland::Projector::D1R1<backend_t>>(Backward::D1R1::id());
      this->mImpl.addOperator<Fft::Worland::Projector::DivR1D1R1_Zero<backend_t>>(Backward::Overr1D1R1::id());
      this->mImpl.addOperator<Fft::Worland::Projector::SphLapl<backend_t>>(Backward::Slapl::id());

      // Create integrators
      this->mImpl.addOperator<Fft::Worland::Integrator::P<backend_t>>(Forward::P::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::R1_Zero<backend_t>>(Forward::Pol::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::DivR1_Zero<backend_t>>(Forward::Q::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::DivR1D1R1_Zero<backend_t>>(Forward::S::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::P_Zero<backend_t>>(Forward::T::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I2<backend_t>>(Forward::I2P::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I2DivR1_Zero<backend_t>>(Forward::I2Q::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I2DivR1D1R1_Zero<backend_t>>(Forward::I2S::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I2_Zero<backend_t>>(Forward::I2T::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I4DivR1_Zero<backend_t>>(Forward::I4Q::id());
      this->mImpl.addOperator<Fft::Worland::Integrator::I4DivR1D1R1_Zero<backend_t>>(Forward::I4S::id());

      // Create reductors
      this->mImpl.addOperator<Fft::Worland::Reductor::EnergySLaplR2<backend_t>>(Reductor::EnergySlaplR2::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::EnergyD1R1<backend_t>>(Reductor::EnergyD1R1::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::EnergyR2<backend_t>>(Reductor::EnergyR2::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::Energy<backend_t>>(Reductor::Energy::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::PowerSLaplR2<backend_t>>(Reductor::PowerSlaplR2::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::PowerD1R1<backend_t>>(Reductor::PowerD1R1::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::PowerR2<backend_t>>(Reductor::PowerR2::id());
      this->mImpl.addOperator<Fft::Worland::Reductor::Power<backend_t>>(Reductor::Power::id());
   }

   void SphereFftWorlandTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereFftWorlandTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereFftWorlandTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void SphereFftWorlandTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFftWorlandTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat SphereFftWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void SphereFftWorlandTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRASPHEREWORLAND, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }
}
}
