/**
 * file SphereWorlandTransform.cpp
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
#include "QuICC/Transform/SphereWorlandTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Poly/Worland/Projector/P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/SphLapl.hpp"

#include "QuICC/Transform/Poly/Worland/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/R1.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2DivR1D1R1_Zero.hpp"

#include "QuICC/Transform/Poly/Worland/Reductor/EnergySLAPLR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Energy.hpp"

#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLAPLR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Power.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/R_1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/R_1D1R1.hpp"
#include "QuICC/Transform/Backward/SLapl.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Q4.hpp"
#include "QuICC/Transform/Forward/S4.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Forward/Q2.hpp"
#include "QuICC/Transform/Forward/S2.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/EnergySLAPLR2.hpp"

#include "QuICC/Transform/Reductor/Power.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Reductor/PowerSLAPLR2.hpp"

#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   SphereWorlandTransform::SphereWorlandTransform()
   {
   }

   SphereWorlandTransform::~SphereWorlandTransform()
   {
   }

   void SphereWorlandTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void SphereWorlandTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array SphereWorlandTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void SphereWorlandTransform::init(SphereWorlandTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the ID to operator mapping
      this->initOperators();
   }

   void SphereWorlandTransform::initOperators()
   {
      // Create projectors
      this->mImpl.addOperator<Poly::Worland::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Poly::Worland::Projector::DivR1_Zero>(Backward::R_1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1R1>(Backward::D1R1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::DivR1D1R1_Zero>(Backward::R_1D1R1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::SphLapl>(Backward::SLapl::id());

      // Create integrators
      this->mImpl.addOperator<Poly::Worland::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::R1>(Forward::R1::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I4DivR1_Zero>(Forward::Q4::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I4DivR1D1R1_Zero>(Forward::S4::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I2_Zero>(Forward::T::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I2DivR1_Zero>(Forward::Q2::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I2DivR1D1R1_Zero>(Forward::S2::id());

      // Create reductors
      this->mImpl.addOperator<Poly::Worland::Reductor::EnergySLAPLR2>(Reductor::EnergySLAPLR2::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::EnergyD1R1>(Reductor::EnergyD1R1::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::EnergyR2>(Reductor::EnergyR2::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::Energy>(Reductor::Energy::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::PowerSLAPLR2>(Reductor::PowerSLAPLR2::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::PowerD1R1>(Reductor::PowerD1R1::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::PowerR2>(Reductor::PowerR2::id());
      this->mImpl.addOperator<Poly::Worland::Reductor::Power>(Reductor::Power::id());
   }

   void SphereWorlandTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereWorlandTransform::forward");
      this->mImpl.transform(rOut, in, id);
   }

   void SphereWorlandTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereWorlandTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void SphereWorlandTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereWorlandTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat SphereWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void SphereWorlandTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRASPHEREWORLAND, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }
}
}
