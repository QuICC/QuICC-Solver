/**
 * @file ALegendreTransform.cpp
 * @brief Source of the implementation of the associated Legendre transform
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/ALegendreTransform.hpp"

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
#include "QuICC/Transform/Forward/OverlaplhD1.hpp"
#include "QuICC/Transform/Forward/Oversin.hpp"
#include "QuICC/Transform/Forward/LaplhOversin.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversin.hpp"
#include "QuICC/Transform/Forward/Laplh.hpp"
#include "QuICC/Transform/Forward/Laplh2.hpp"
#include "QuICC/Transform/Forward/Overlaplh.hpp"
#include "QuICC/Transform/Forward/OversinDphi.hpp"
#include "QuICC/Transform/Forward/LaplhOversinDphi.hpp"
#include "QuICC/Transform/Forward/OverlaplhOversinDphi.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1Laplh.hpp"
#include "QuICC/Transform/Backward/Oversin.hpp"
#include "QuICC/Transform/Backward/OversinDphi.hpp"
#include "QuICC/Transform/Backward/OversinLaplh.hpp"
#include "QuICC/Transform/Backward/OversinLaplhDphi.hpp"
#include "QuICC/Transform/Backward/OversinD1Sin.hpp"

namespace QuICC {

namespace Transform {

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
      using namespace Poly::ALegendre;
      // Reserve storage for the projectors, 1/sin projectors and derivative
#ifdef QUICC_USE_KOKKOS_CUDA
      this->mImpl.addOperator<Poly::ALegendre::Projector::P<kokkos_t>>(Backward::P::id());
      using backend_t = kokkos_t;
#elif QUICC_USE_KOKKOS
      this->mImpl.addOperator<Poly::ALegendre::Projector::P<kokkos_t>>(Backward::P::id());
      using backend_t = base_t;
#else
      this->mImpl.addOperator<Poly::ALegendre::Projector::P<base_t>>(Backward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1D1S1<base_t>>(Backward::OversinD1Sin::id());
      using backend_t = base_t;	
#endif
   this->mImpl.addOperator<Poly::ALegendre::Projector::Ll<backend_t>>(Backward::Laplh::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::D1<backend_t>>(Backward::D1::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::LlD1<backend_t>>(Backward::D1Laplh::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1<backend_t>>(Backward::Oversin::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1Dp<backend_t>>(Backward::OversinDphi::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1<backend_t>>(Backward::OversinLaplh::id());
   this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1Dp<backend_t>>(Backward::OversinLaplhDphi::id());

      // Reserve storage for the weighted projectors, 1/sin projectors and derivative
#ifdef QUICC_USE_KOKKOS_CUDA
      this->mImpl.addOperator<Poly::ALegendre::Integrator::P<kokkos_t>>(Forward::P::id());
      using backend_t = kokkos_t;
#elif QUICC_USE_KOKKOS
      this->mImpl.addOperator<Poly::ALegendre::Integrator::P<kokkos_t>>(Forward::P::id());
      using backend_t = base_t;
#else
      this->mImpl.addOperator<Poly::ALegendre::Integrator::P<base_t>>(Forward::P::id());
      using backend_t = base_t;
#endif
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLl<backend_t>>(Forward::Overlaplh::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll<backend_t>>(Forward::Laplh::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll2<backend_t>>(Forward::Laplh2::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1<backend_t>>(Forward::Oversin::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1<backend_t>>(Forward::OverlaplhOversin::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1<backend_t>>(Forward::LaplhOversin::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::D1<backend_t>>(Forward::D1::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlD1<backend_t>>(Forward::OverlaplhD1::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::LlD1<backend_t>>(Forward::LaplhD1::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1Dp<backend_t>>(Forward::OversinDphi::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1Dp<backend_t>>(Forward::OverlaplhOversinDphi::id());
   this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1Dp<backend_t>>(Forward::LaplhOversinDphi::id());

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
