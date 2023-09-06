/** 
 * @file TTT.cpp
 * @brief Source of TTT of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme implementation
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/TTT.hpp"
#include "QuICC/SpatialScheme/3D/TTTBuilder.hpp"
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/FullyCoupled.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string TTT::sTag = "TTT";

   const std::string TTT::sFormatted = "TTT";

   const std::size_t TTT::sId = Hasher::makeId(TTT::sTag);

   TTT::TTT(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, TTT::sId, TTT::sTag, TTT::sFormatted)
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
      this->enable(Feature::RegularSpectrum);
      this->enable(Feature::SpectralMatrix3D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::TransformSpectralOrdering123);
      this->enable(Feature::RealSpectrum);
   }

   std::shared_ptr<IBuilder> TTT::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<TTTBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> TTT::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
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

   std::shared_ptr<Parallel::IIndexConv> TTT::createIndexConv(const Dimensions::Transform::Id id) const
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

   std::shared_ptr<Equations::Tools::ICoupling> TTT::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SINGLE)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::FullyCoupled>();
      return spCoupling;
   }

   TTT::VariantTransformDataPointer TTT::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TTT::VariantTransformDataPointer TTT::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<TTT::RealTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   TTT::ScalarVariable TTT::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      TTT::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   TTT::VectorVariable TTT::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      TTT::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

} // SpatialScheme
} // QuICC
