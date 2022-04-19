/** 
 * @file ALegendreTransform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/ALegendreTransform.hpp"

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/P.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/Ll.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/D1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/LlD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/LlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/LlDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/DivS1D1S1.hpp"
                                         
#include "QuICC/Transform/Poly/ALegendre/Integrator/P.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Ll.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/Ll2.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLl.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/LlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLlDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/LlDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/D1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/DivLlD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/LlD1.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/D1.hpp"
#include "QuICC/Transform/Forward/LaplhD1.hpp"
#include "QuICC/Transform/Forward/Laplh_1D1.hpp"
#include "QuICC/Transform/Forward/Sin_1.hpp"
#include "QuICC/Transform/Forward/LaplhSin_1.hpp"
#include "QuICC/Transform/Forward/Laplh_1Sin_1.hpp"
#include "QuICC/Transform/Forward/Laplh.hpp"
#include "QuICC/Transform/Forward/Laplh2.hpp"
#include "QuICC/Transform/Forward/Laplh_1.hpp"
#include "QuICC/Transform/Forward/Sin_1Dphi.hpp"
#include "QuICC/Transform/Forward/LaplhSin_1Dphi.hpp"
#include "QuICC/Transform/Forward/Laplh_1Sin_1Dphi.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1Laplh.hpp"
#include "QuICC/Transform/Backward/Sin_1.hpp"
#include "QuICC/Transform/Backward/Sin_1Dphi.hpp"
#include "QuICC/Transform/Backward/Sin_1Laplh.hpp"
#include "QuICC/Transform/Backward/Sin_1LaplhDphi.hpp"
#include "QuICC/Transform/Backward/Sin_1D1Sin.hpp"
                                         
namespace QuICC {

namespace Transform {

   ALegendreTransform::ALegendreTransform()
   {
   }

   ALegendreTransform::~ALegendreTransform()
   {
   }

   void ALegendreTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void ALegendreTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array ALegendreTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void ALegendreTransform::init(ALegendreTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the quadrature grid and weights and operators
      this->initOperators();
   }

   void ALegendreTransform::initOperators()
   {
      // Reserve storage for the projectors, 1/sin projectors and derivative
      this->mImpl.addOperator<Poly::ALegendre::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::Ll>(Backward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlD1>(Backward::D1Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1>(Backward::Sin_1::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1Dp>(Backward::Sin_1Dphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1>(Backward::Sin_1Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1Dp>(Backward::Sin_1LaplhDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1D1S1>(Backward::Sin_1D1Sin::id());

      // Reserve storage for the weighted projectors, 1/sin projectors and derivative
      this->mImpl.addOperator<Poly::ALegendre::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLl>(Forward::Laplh_1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll>(Forward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll2>(Forward::Laplh2::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1>(Forward::Sin_1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1>(Forward::Laplh_1Sin_1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1>(Forward::LaplhSin_1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::D1>(Forward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlD1>(Forward::Laplh_1D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlD1>(Forward::LaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1Dp>(Forward::Sin_1Dphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1Dp>(Forward::Laplh_1Sin_1Dphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1Dp>(Forward::LaplhSin_1Dphi::id());
   }

   void ALegendreTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void ALegendreTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void ALegendreTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::reduce(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ALegendreTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat ALegendreTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void ALegendreTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRAALEGENDRE, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
