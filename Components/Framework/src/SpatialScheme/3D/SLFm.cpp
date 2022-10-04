/** 
 * @file SLFm.cpp
 * @brief Source of the spherical Chebyshev(FFT) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with m spectral ordering
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/SLFm.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/3D/SLFmBuilder.hpp"
#include "QuICC/Transform/ShellChebyshevTransform.hpp"
#include "QuICC/Transform/ALegendreTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/SHmIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/Equations/Tools/SHm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string SLFm::sTag = "SLFm";

   const std::string SLFm::sFormatted = "SLFm";

   const std::size_t SLFm::sId = registerId<SLFm>(SLFm::sTag);

   SLFm::SLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, SLFm::sId, SLFm::sTag, SLFm::sFormatted)
   {
      // pHYSIcal component aliases
      this->mPhys.add(FieldComponents::Physical::R);
      this->mPhys.add(FieldComponents::Physical::THETA);
      this->mPhys.add(FieldComponents::Physical::PHI);

      // Spectral component aliases
      if(formulation == VectorFormulation::TORPOL)
      {
         this->mSpec.add(FieldComponents::Spectral::TOR);
         this->mSpec.add(FieldComponents::Spectral::POL);
         this->mSpec.add(FieldComponents::Spectral::NOTUSED);
      } else if(formulation == VectorFormulation::QST)
      {
         this->mSpec.add(FieldComponents::Spectral::Q);
         this->mSpec.add(FieldComponents::Spectral::S);
         this->mSpec.add(FieldComponents::Spectral::T);
      }

      // Enable basic scheme features
      this->enable(Feature::ShellGeometry);
      this->enable(Feature::FourierIndex3);
      this->enable(Feature::SpectralMatrix2D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::TransformSpectralOrdering123);
      this->enable(Feature::ComplexSpectrum);
   }

   SLFm::~SLFm()
   {
   }

   std::shared_ptr<IBuilder> SLFm::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<SLFmBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> SLFm::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spTT = std::make_shared<Transform::ShellChebyshevTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::ShellChebyshevTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spTT->init(spST);
            spTransform = spTT;
            break;
         }
         case Dimensions::Transform::TRA2D:
         {
            auto spAT = std::make_shared<Transform::ALegendreTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::ALegendreTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spAT->init(spST);
            spTransform = spAT;
            break;
         }
         case Dimensions::Transform::TRA3D:
         {
            auto spFT = std::make_shared<Transform::MixedFourierTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::MixedFourierTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spFT->init(spST);
            spTransform = spFT;
            break;
         }
         default:
            throw std::logic_error("Tried to initialize too many transforms");
      }

      return spTransform;
   }

   std::shared_ptr<Parallel::IIndexConv> SLFm::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            spConv = std::make_shared<Parallel::PassthroughIndexConv>();
            break;
         case Dimensions::Transform::TRA2D:
            spConv = std::make_shared<Parallel::SHmIndexConv>();
            break;
         case Dimensions::Transform::TRA3D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         default:
            throw std::logic_error("Tried to initialize an impossible index converter");
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> SLFm::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::SHm>();
      return spCoupling;
   }

   SLFm::VariantTransformDataPointer SLFm::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<SLFm::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   SLFm::VariantTransformDataPointer SLFm::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<SLFm::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   SLFm::ScalarVariable SLFm::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      SLFm::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   SLFm::VectorVariable SLFm::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      SLFm::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
