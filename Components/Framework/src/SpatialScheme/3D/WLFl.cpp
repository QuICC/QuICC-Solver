/** 
 * @file WLFl.cpp
 * @brief Source of the ID for a sphere Worland(poly) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/3D/WLFlBuilder.hpp"
#include "QuICC/SpatialScheme/3D/WLFlfftBuilder.hpp"
#include "QuICC/SpatialScheme/3D/WLFlTriangularBuilder.hpp"
#include "QuICC/Transform/SphereWorlandTransform.hpp"
#include "QuICC/Transform/SphereFftWorlandTransform.hpp"
#include "QuICC/Transform/ALegendreTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/SHlIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/Equations/Tools/SHl.hpp"
#include "QuICC/Equations/Tools/SHlm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string WLFl::sTag = "WLFl";

   const std::string WLFl::sFormatted = "WLFl";

   const std::size_t WLFl::sId = registerId<WLFl>(WLFl::sTag);

   WLFl::WLFl(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, WLFl::sId, WLFl::sTag, WLFl::sFormatted)
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
      }
      else if(formulation == VectorFormulation::QST)
      {
         this->mSpec.add(FieldComponents::Spectral::Q);
         this->mSpec.add(FieldComponents::Spectral::S);
         this->mSpec.add(FieldComponents::Spectral::T);
      }

      // Enable basic scheme features
      this->enable(Feature::SphereGeometry);
      this->enable(Feature::FourierIndex2);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering132);
      this->enable(Feature::TransformSpectralOrdering132);
      this->enable(Feature::ComplexSpectrum);
   }

   std::shared_ptr<IBuilder> WLFl::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      std::shared_ptr<IBuilder>  spIB;
      if(this->mImplType(0) == 1)
      {
         auto spBuilder = makeBuilder<WLFlfftBuilder>(dim, this->purpose(), needInterpretation);
         spIB = spBuilder;
      }
      else if(this->mImplType(0) == 2)
      {
         auto spBuilder = makeBuilder<WLFlTriangularBuilder>(dim, this->purpose(), needInterpretation);
         spIB = spBuilder;
      } else
      {
         auto spBuilder = makeBuilder<WLFlBuilder>(dim, this->purpose(), needInterpretation);
         spIB = spBuilder;
      }

      return spIB;
   }

   std::shared_ptr<Transform::ITransform> WLFl::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            if(this->mImplType(0) == 1)
            {
               auto spWT = std::make_shared<Transform::SphereFftWorlandTransform>();
               auto spST = std::dynamic_pointer_cast<Transform::SphereFftWorlandTransform::SetupType>(spSetup);
               if(!spST)
               {
                  throw std::logic_error("Incompatible transform setup given");
               }
               spWT->init(spST);
               spTransform = spWT;
            } else
            {
               auto spWT = std::make_shared<Transform::SphereWorlandTransform>();
               auto spST = std::dynamic_pointer_cast<Transform::SphereWorlandTransform::SetupType>(spSetup);
               if(!spST)
               {
                  throw std::logic_error("Incompatible transform setup given");
               }
               spWT->init(spST);
               spTransform = spWT;
            }
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

   std::shared_ptr<Parallel::IIndexConv> WLFl::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            spConv = std::make_shared<Parallel::PassthroughIndexConv>();
            break;
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

   std::shared_ptr<Equations::Tools::ICoupling> WLFl::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType == Equations::CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         auto spCoupling = std::make_shared<Equations::Tools::SHl>();
         return spCoupling;
      }
      else if(indexType == Equations::CouplingIndexType::MODE)
      {
         auto spCoupling = std::make_shared<Equations::Tools::SHlm>();
         return spCoupling;
      }
      else
      {
         throw std::logic_error("Unknown coupling tools");
      }
   }

   WLFl::VariantTransformDataPointer WLFl::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<WLFl::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   WLFl::VariantTransformDataPointer WLFl::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<WLFl::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   WLFl::ScalarVariable WLFl::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      WLFl::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   WLFl::VectorVariable WLFl::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      WLFl::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

} // SpatialScheme
} // QuICC
