/** 
 * @file TFT.cpp
 * @brief Source of ID of the Chebyshev(FFT) + Fourier + Chebyshev(FFT) scheme implementation
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/TFT.hpp"
#include "QuICC/SpatialScheme/3D/TFTBuilder.hpp"
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/Regular1D.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string TFT::sTag = "TFT";

   const std::string TFT::sFormatted = "TFT";

   const std::size_t TFT::sId = Hasher::makeId(TFT::sTag);

   TFT::TFT(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, TFT::sId, TFT::sTag, TFT::sFormatted)
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
      this->enable(Feature::FourierIndex3);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix2D);
      this->enable(Feature::SpectralOrdering132);
      this->enable(Feature::TransformSpectralOrdering132);
      this->enable(Feature::ComplexSpectrum);
   }

   std::shared_ptr<IBuilder> TFT::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<TFTBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> TFT::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
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
         case Dimensions::Transform::TRA3D:
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
         default:
            throw std::logic_error("Tried to initialize too many transforms");
      }

      return spTransform;
   }

   std::shared_ptr<Parallel::IIndexConv> TFT::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA2D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         case Dimensions::Transform::TRA3D:
            spConv = std::make_shared<Parallel::NoIndexConv>();
            break;
         default:
            throw std::logic_error("Tried to initialize an impossible index converter");
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> TFT::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::Regular1D>();
      return spCoupling;
   }

   TFT::VariantTransformDataPointer TFT::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TFT::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TFT::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TFT::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TFT::VariantTransformDataPointer TFT::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TFT::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TFT::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TFT::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TFT::ScalarVariable TFT::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      TFT::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   TFT::VectorVariable TFT::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      TFT::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

} // SpatialScheme
} // QuICC
