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
#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Projector/PP.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PLl.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PLlD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PLlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Projector/PLlDivS1Dp.hpp"
/* #include "QuICC/Transform/Poly/ALegendre/Projector/PDivS1D1S1.hpp" */
#endif

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

#ifdef QUICC_USE_KOKKOS
#include "QuICC/Transform/Poly/ALegendre/Integrator/PP.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLl.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLl2.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLl.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLlDivS1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLlDivS1Dp.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PDivLlD1.hpp"
#include "QuICC/Transform/Poly/ALegendre/Integrator/PLlD1.hpp"
#endif

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
#ifdef QUICC_USE_KOKKOS_CUDA
      this->mImpl.addOperator<Poly::ALegendre::Projector::PP<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PLl<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PD1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PLlD1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::D1Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PDivS1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PDivS1Dp<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PLlDivS1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::OversinLaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::PLlDivS1Dp<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::OversinLaplhDphi::id());

      /* this->mImpl.addOperator<Poly::ALegendre::Projector::PDivS1D1S1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Backward::OversinD1Sin::id()); */
#elif QUICC_USE_KOKKOS
      this->mImpl.addOperator<Poly::ALegendre::Projector::PP<Poly::ALegendre::PIALegendreOperatorTypes>>(Backward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::Ll>(Backward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlD1>(Backward::D1Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1>(Backward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1Dp>(Backward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1>(Backward::OversinLaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1Dp>(Backward::OversinLaplhDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1D1S1>(Backward::OversinD1Sin::id());
#else
      this->mImpl.addOperator<Poly::ALegendre::Projector::P<>>(Backward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::Ll>(Backward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlD1>(Backward::D1Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1>(Backward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1Dp>(Backward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1>(Backward::OversinLaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::LlDivS1Dp>(Backward::OversinLaplhDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Projector::DivS1D1S1>(Backward::OversinD1Sin::id());
#endif


      // Reserve storage for the weighted projectors, 1/sin projectors and derivative
#ifdef QUICC_USE_KOKKOS_CUDA
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PP<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivLl<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::Overlaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PLl<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PLl2<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::Laplh2::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivS1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivLlDivS1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::OverlaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PLlDivS1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::LaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PD1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivLlD1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::OverlaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PLlD1<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::LaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivS1Dp<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PDivLlDivS1Dp<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::OverlaplhOversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PLlDivS1Dp<Poly::ALegendre::CudaIALegendreOperatorTypes>>(Forward::LaplhOversinDphi::id());
#elif QUICC_USE_KOKKOS
      this->mImpl.addOperator<Poly::ALegendre::Integrator::PP<Poly::ALegendre::PIALegendreOperatorTypes>>(Forward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLl>(Forward::Overlaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll>(Forward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll2>(Forward::Laplh2::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1>(Forward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1>(Forward::OverlaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1>(Forward::LaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::D1>(Forward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlD1>(Forward::OverlaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlD1>(Forward::LaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1Dp>(Forward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1Dp>(Forward::OverlaplhOversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1Dp>(Forward::LaplhOversinDphi::id());
#else
      this->mImpl.addOperator<Poly::ALegendre::Integrator::P<>>(Forward::P::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLl>(Forward::Overlaplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll>(Forward::Laplh::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::Ll2>(Forward::Laplh2::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1>(Forward::Oversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1>(Forward::OverlaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1>(Forward::LaplhOversin::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::D1>(Forward::D1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlD1>(Forward::OverlaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlD1>(Forward::LaplhD1::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivS1Dp>(Forward::OversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::DivLlDivS1Dp>(Forward::OverlaplhOversinDphi::id());
      this->mImpl.addOperator<Poly::ALegendre::Integrator::LlDivS1Dp>(Forward::LaplhOversinDphi::id());
#endif

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
