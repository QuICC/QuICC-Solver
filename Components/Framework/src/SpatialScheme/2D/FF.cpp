/**
 * @file FF.cpp
 * @brief Source of the Fourier + Fourier scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/2D/FF.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/2D/FFBuilder.hpp"
#include "QuICC/Transform/ComplexFourierTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string FF::sTag = "FF";

   const std::string FF::sFormatted = "FF";

   const std::size_t FF::sId = registerId<FF>(FF::sTag);

   FF::FF(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 2, FF::sId, FF::sTag, FF::sFormatted)
   {
      // Physical component aliases
      this->mPhys.add(FieldComponents::Physical::X);
      this->mPhys.add(FieldComponents::Physical::Y);
      this->mPhys.add(FieldComponents::Physical::NOTUSED);

      // Spectral component aliases
      this->mSpec.add(FieldComponents::Spectral::X);
      this->mSpec.add(FieldComponents::Spectral::Y);
      this->mSpec.add(FieldComponents::Spectral::NOTUSED);

      // Enable basic scheme features
      this->enable(Feature::CartesianGeometry);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::ComplexSpectrum);
   }

   FF::~FF()
   {
   }

   std::shared_ptr<IBuilder> FF::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<FFBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> FF::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spTT = std::make_shared<Transform::ComplexFourierTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::ComplexFourierTransform::SetupType>(spSetup);
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

   std::shared_ptr<Parallel::IIndexConv> FF::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA2D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         default:
            throw std::logic_error("Tried to initialize an impossible index converter");
      }

      return spConv;
   }

   FF::VariantTransformDataPointer FF::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<FF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<FF::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   FF::VariantTransformDataPointer FF::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<FF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<FF::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   FF::ScalarVariable FF::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      FF::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   FF::VectorVariable FF::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      FF::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
