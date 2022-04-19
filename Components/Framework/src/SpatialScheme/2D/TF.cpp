/**
 * @file TF.cpp
 * @brief Source of the Chebyshev(FFT) + Fourier scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/2D/TF.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Coordinator.hpp"
#include "QuICC/SpatialScheme/2D/TFBuilder.hpp"
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/Regular1D.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string TF::sTag = "TF";

   const std::string TF::sFormatted = "TF";

   const std::size_t TF::sId = registerId<TF>(TF::sTag);

   TF::TF(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 2, TF::sId, TF::sTag, TF::sFormatted)
   {
      // Physical component aliases
      this->mPhys.add(FieldComponents::Physical::Z);
      this->mPhys.add(FieldComponents::Physical::X);
      this->mPhys.add(FieldComponents::Physical::NOTUSED);

      // Spectral component aliases
      this->mSpec.add(FieldComponents::Spectral::Z);
      this->mSpec.add(FieldComponents::Spectral::X);
      this->mSpec.add(FieldComponents::Spectral::NOTUSED);

      // Enable basic scheme features
      this->enable(Feature::CartesianGeometry);
      this->enable(Feature::FourierIndex2);
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::ComplexSpectrum);
   }

   TF::~TF()
   {
   }

   std::shared_ptr<IBuilder> TF::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<TFBuilder>(dim, this->purpose(), needInterpretation);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> TF::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
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
         default:
            throw std::logic_error("Tried to initialize too many transforms");
      }

      return spTransform;
   }

   std::shared_ptr<Parallel::IIndexConv> TF::createIndexConv(const Dimensions::Transform::Id id) const
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

   std::shared_ptr<Equations::Tools::ICoupling> TF::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::Regular1D>();
      return spCoupling;
   }

   TF::VariantTransformDataPointer TF::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TF::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TF::VariantTransformDataPointer TF::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TF::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TF::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TF::ScalarVariable TF::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      TF::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   TF::VectorVariable TF::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      TF::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

}
}
