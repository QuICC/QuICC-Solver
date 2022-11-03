/** 
 * @file CylinderWorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform in a cylinder
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/CylinderWorlandTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Poly/Worland/Projector/P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1_P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/CylLaplh.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/CylLaplh_DivR1D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1CylLaplh_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1CylLaplh.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1CylLaplh_D1DivR1D1R1.hpp"

#include "QuICC/Transform/Poly/Worland/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1D1R1_I2.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6DivR1D1R1_I4.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6CylLaplh_I4D1R1.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I4Overr1Pm.hpp"
#include "QuICC/Transform/Forward/I4Overr1D1R1ZI2.hpp"
#include "QuICC/Transform/Forward/I6Overr1Pm.hpp"
#include "QuICC/Transform/Forward/I6Overr1D1R1ZI4.hpp"
#include "QuICC/Transform/Forward/I6LaplhZI4D1R1.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1Pm.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1ZP.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"
#include "QuICC/Transform/Backward/LaplhZOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/Overr1LaplhPm.hpp"
#include "QuICC/Transform/Backward/D1Laplh.hpp"
#include "QuICC/Transform/Backward/D1LaplhZD1Overr1D1R1.hpp"

namespace QuICC {

namespace Transform {

   CylinderWorlandTransform::CylinderWorlandTransform()
   {
   }

   CylinderWorlandTransform::~CylinderWorlandTransform()
   {
   }

   void CylinderWorlandTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void CylinderWorlandTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array CylinderWorlandTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void CylinderWorlandTransform::init(CylinderWorlandTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the quadrature grid and weights and operators
      this->initOperators();
   }

   void CylinderWorlandTransform::initOperators()
   {
      // Reserve storage for the projectors
      this->mImpl.addOperator<Poly::Worland::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Poly::Worland::Projector::DivR1_Zero>(Backward::Overr1Pm::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1_P>(Backward::D1ZP::id());
      this->mImpl.addOperator<Poly::Worland::Projector::DivR1D1R1>(Backward::Overr1D1R1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::CylLaplh>(Backward::Laplh::id());
      this->mImpl.addOperator<Poly::Worland::Projector::CylLaplh_DivR1D1R1>(Backward::LaplhZOverr1D1R1::id());
      this->mImpl.addOperator<Poly::Worland::Projector::DivR1CylLaplh_Zero>(Backward::Overr1LaplhPm::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1CylLaplh>(Backward::D1Laplh::id());
      this->mImpl.addOperator<Poly::Worland::Projector::D1CylLaplh_D1DivR1D1R1>(Backward::D1LaplhZD1Overr1D1R1::id());

      // Reserve storage for the weighted projectors 
      this->mImpl.addOperator<Poly::Worland::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I4DivR1_Zero>(Forward::I4Overr1Pm::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I4DivR1D1R1_I2>(Forward::I4Overr1D1R1ZI2::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I6DivR1_Zero>(Forward::I6Overr1Pm::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I6DivR1D1R1_I4>(Forward::I6Overr1D1R1ZI4::id());
      this->mImpl.addOperator<Poly::Worland::Integrator::I6CylLaplh_I4D1R1>(Forward::I6LaplhZI4D1R1::id());
   }

   void CylinderWorlandTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CylinderWorlandTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void CylinderWorlandTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::reduce(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CylinderWorlandTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat CylinderWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void CylinderWorlandTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRACYLINDERWORLAND, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
