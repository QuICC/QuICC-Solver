/** 
 * @file TFF.cpp
 * @brief Source of ID of the Chebyshev(FFT) + Fourier + Fourier scheme implementation
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/TFF.hpp"
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/SpatialScheme/3D/TFFBuilder.hpp"
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Transform/ComplexFourierTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/PMIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/Equations/Tools/Regular2D.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string TFF::sTag = "TFF";

   const std::string TFF::sFormatted = "TFF";

   const std::size_t TFF::sId = Hasher::makeId(TFF::sTag);

   TFF::TFF(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, TFF::sId, TFF::sTag, TFF::sFormatted)
   {
      // Physical component aliases
      this->mPhys.add(FieldComponents::Physical::Z);
      this->mPhys.add(FieldComponents::Physical::X);
      this->mPhys.add(FieldComponents::Physical::Y);

      // Spectral component aliases
      if(formulation == VectorFormulation::TORPOL)
      {
         this->mSpec.add(FieldComponents::Spectral::TOR);
         this->mSpec.add(FieldComponents::Spectral::POL);
         this->mSpec.add(FieldComponents::Spectral::NOTUSED);
      }
      else
      {
         this->mSpec.add(FieldComponents::Spectral::Z);
         this->mSpec.add(FieldComponents::Spectral::X);
         this->mSpec.add(FieldComponents::Spectral::Y);
      }

      // Enable basic scheme features
      this->enable(Feature::CartesianGeometry);
      this->enable(Feature::FourierIndex23);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering132);
      this->enable(Feature::TransformSpectralOrdering132);
      this->enable(Feature::ComplexSpectrum);
   }

   std::shared_ptr<IBuilder> TFF::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<TFFBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> TFF::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spTT = std::make_shared<Transform::CartesianChebyshevTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::CartesianChebyshevTransform::SetupType>(spSetup);
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

   std::shared_ptr<Parallel::IIndexConv> TFF::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            spConv = std::make_shared<Parallel::PassthroughIndexConv>();
            break;
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

   std::shared_ptr<Equations::Tools::ICoupling> TFF::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::MODE)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::Regular2D>();
      return spCoupling;
   }

   TFF::VariantTransformDataPointer TFF::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TFF::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TFF::VariantTransformDataPointer TFF::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<TFF::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TFF::ScalarVariable TFF::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      TFF::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   TFF::VectorVariable TFF::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      TFF::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

} // SpatialScheme
} // QuICC
