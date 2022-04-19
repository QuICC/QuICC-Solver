/** 
 * @file SLFl.cpp
 * @brief Source of ID of the spherical Chebyshev(FFT) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with l spectral ordering
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/SLFl.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/3D/SLFlBuilder.hpp"
#include "QuICC/Transform/ShellChebyshevTransform.hpp"
#include "QuICC/Transform/ALegendreTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/SHlIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/SHl.hpp"
#include "QuICC/Equations/Tools/SHlm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string SLFl::sTag = "SLFl";

   const std::string SLFl::sFormatted = "SLFl";

   const std::size_t SLFl::sId = registerId<SLFl>(SLFl::sTag);

   SLFl::SLFl(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, SLFl::sId, SLFl::sTag, SLFl::sFormatted)
   {
      // Physical component aliases
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
      this->enable(Feature::FourierIndex2);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering132);
      this->enable(Feature::ComplexSpectrum);
   }

   SLFl::~SLFl()
   {
   }

   std::shared_ptr<IBuilder> SLFl::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<SLFlBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> SLFl::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
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

   std::shared_ptr<Parallel::IIndexConv> SLFl::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA2D:
            spConv = std::make_shared<Parallel::SHlIndexConv>();
            break;
         case Dimensions::Transform::TRA3D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         default:
            throw std::logic_error("Tried to initialize an impossible index converter");
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> SLFl::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType == Equations::CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         auto spCoupling = std::make_shared<Equations::Tools::SHl>();
         return spCoupling;
      } else if(indexType == Equations::CouplingIndexType::MODE)
      {
         auto spCoupling = std::make_shared<Equations::Tools::SHlm>();
         return spCoupling;
      } else
      {
         throw std::logic_error("Unknown coupling tools");
      }
   }

   SLFl::VariantTransformDataPointer SLFl::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<SLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<SLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<SLFl::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   SLFl::VariantTransformDataPointer SLFl::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<SLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<SLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<SLFl::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   SLFl::ScalarVariable SLFl::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      SLFl::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   SLFl::VectorVariable SLFl::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      SLFl::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
