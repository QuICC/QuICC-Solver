/** 
 * @file WLFm.cpp
 * @brief Source of ID of the sphere Worland(poly) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral m ordering
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/WLFm.hpp"
#include "QuICC/SpatialScheme/3D/WLFmBuilder.hpp"
#include "QuICC/Transform/SphereWorlandTransform.hpp"
#include "QuICC/Transform/ALegendreTransform.hpp"
#include "QuICC/Transform/MixedFourierTransform.hpp"
#include "QuICC/Transform/Setup/Default.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/Communicators/Converters/SHm2lIndexConv.hpp"
#include "QuICC/Communicators/Converters/SHlIndexConv.hpp"
#include "QuICC/Communicators/Converters/NoIndexConv.hpp"
#include "QuICC/Equations/Tools/SHm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string WLFm::sTag = "WLFm";

   const std::string WLFm::sFormatted = "WLFm";

   const std::size_t WLFm::sId = Hasher::makeId(WLFm::sTag);

   WLFm::WLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : ISpatialScheme(formulation, purpose, 3, WLFm::sId, WLFm::sTag, WLFm::sFormatted)
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
      this->enable(Feature::FourierIndex3);
      this->enable(Feature::SpectralMatrix2D);
      this->enable(Feature::SpectralOrdering123);
      this->enable(Feature::TransformSpectralOrdering132);
      this->enable(Feature::ComplexSpectrum);
   }

   void WLFm::setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type)
   {
      assert(type.size() == this->mImplType.size());

      // Replace default with effective implementation options for 1D
      assert(type.size() > 0);
      std::size_t dimId = 0;
      const auto& opt1D = type.at(dimId);
      auto& mOpt1D = this->mImplType.at(dimId);
      mOpt1D.clear();
      if(std::find(opt1D.begin(), opt1D.end(), Transform::Setup::Default::id()) != opt1D.end())
      {
         mOpt1D.push_back(Transform::Setup::GaussianQuadrature::id());
         mOpt1D.push_back(Transform::Setup::Uniform::id());
      }
      else
      {
         mOpt1D = opt1D;
      }

      // Replace default with effective implementation options for 2D
      assert(type.size() > 0);
      dimId = 1;
      const auto& opt2D = type.at(dimId);
      auto& mOpt2D = this->mImplType.at(dimId);
      mOpt2D.clear();
      if(std::find(opt2D.begin(), opt2D.end(), Transform::Setup::Default::id()) != opt2D.end())
      {
         mOpt2D.push_back(Transform::Setup::GaussianQuadrature::id());
      }
      else
      {
         mOpt2D = opt2D;
      }

      // Replace default with effective implementation options for 3D
      assert(type.size() > 0);
      dimId = 2;
      const auto& opt3D = type.at(dimId);
      auto& mOpt3D = this->mImplType.at(dimId);
      mOpt3D.clear();
      if(std::find(opt3D.begin(), opt3D.end(), Transform::Setup::Default::id()) != opt3D.end())
      {
         mOpt3D.push_back(Transform::Setup::Fft::id());
      }
      else
      {
         mOpt3D = opt3D;
      }
   }

   std::shared_ptr<IBuilder> WLFm::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<WLFmBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> WLFm::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            auto spWT = std::make_shared<Transform::SphereWorlandTransform>();
            auto spST = std::dynamic_pointer_cast<Transform::SphereWorlandTransform::SetupType>(spSetup);
            if(!spST)
            {
               throw std::logic_error("Incompatible transform setup given");
            }
            spWT->init(spST);
            spTransform = spWT;
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

   std::shared_ptr<Parallel::IIndexConv> WLFm::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            spConv = std::make_shared<Parallel::SHm2lIndexConv>();
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

   std::shared_ptr<Equations::Tools::ICoupling> WLFm::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::SHm>();
      return spCoupling;
   }

   WLFm::VariantTransformDataPointer WLFm::fwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<WLFm::RealTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   WLFm::VariantTransformDataPointer WLFm::bwdPtr(const Dimensions::Transform::Id id) const
   {
      VariantTransformDataPointer v;
      switch(id) {
         case Dimensions::Transform::TRA1D:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA2D:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::TRA3D:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         case Dimensions::Transform::SPECTRAL:
            v = std::forward<WLFm::ComplexTransformDataType *>(0);
            break;
         default:
            throw std::logic_error("Requested forward pointer for unknown dimension");
      }

      return v;
   }

   WLFm::ScalarVariable WLFm::createSVar(std::shared_ptr<Resolution> spRes) const
   {
      WLFm::ScalarVariable p = std::make_shared<Framework::Selector::ScalarVariable<MHDComplex> >(spRes);

      return p;
   }

   WLFm::VectorVariable WLFm::createVVar(std::shared_ptr<Resolution> spRes) const
   {
      WLFm::VectorVariable p = std::make_shared<Framework::Selector::VectorVariable<MHDComplex> >(spRes);

      return p;
   }

} // SpatialScheme
} // QuICC
