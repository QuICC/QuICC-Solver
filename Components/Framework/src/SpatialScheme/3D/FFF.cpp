/** 
 * @file FFF.cpp
 * @brief Source of ID of the Fourier + Fourier + Fourier scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/FFF.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/3D/FFFBuilder.hpp"
#include "QuICC/Transform/ComplexFourierTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/PMIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/Regular3D.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string FFF::sTag = "FFF";

   const std::string FFF::sFormatted = "FFF";

   const std::size_t FFF::sId = registerId<FFF>(FFF::sTag);

   FFF::FFF(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, FFF::sId, FFF::sTag, FFF::sFormatted)
   {
      // Physical component aliases
      this->mPhys.add(FieldComponents::Physical::X);
      this->mPhys.add(FieldComponents::Physical::Y);
      this->mPhys.add(FieldComponents::Physical::Z);

      // Spectral component aliases
      this->mSpec.add(FieldComponents::Spectral::X);
      this->mSpec.add(FieldComponents::Spectral::Y);
      this->mSpec.add(FieldComponents::Spectral::Z);

      // Enable basic scheme features
      this->enable(Feature::CartesianGeometry);
      this->enable(Feature::FourierIndex123);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::TransformSpectralOrdering123);
      this->enable(Feature::ComplexSpectrum);
   }

   FFF::~FFF()
   {
   }

   std::shared_ptr<IBuilder> FFF::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<FFFBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> FFF::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spFT = std::make_shared<Transform::ComplexFourierTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::ComplexFourierTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spFT->init(spST);
            spTransform = spFT;
            break;
         }
         case Dimensions::Transform::TRA2D:
         {
            auto spFT = std::make_shared<Transform::ComplexFourierTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::ComplexFourierTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spFT->init(spST);
            spTransform = spFT;
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

   std::shared_ptr<Parallel::IIndexConv> FFF::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA2D:
            spConv = std::make_shared<Parallel::PMIndexConv>();
            break;
         case Dimensions::Transform::TRA3D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         default:
            throw std::logic_error("Tried to initialize an impossible index converter");
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> FFF::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SINGLE)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::Regular3D>();
      return spCoupling;
   }

   FFF::VariantTransformDataPointer FFF::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<FFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<FFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<FFF::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   FFF::VariantTransformDataPointer FFF::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<FFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<FFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<FFF::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   FFF::ScalarVariable FFF::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      FFF::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   FFF::VectorVariable FFF::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      FFF::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
