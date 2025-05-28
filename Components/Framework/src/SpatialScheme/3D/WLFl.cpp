/** 
 * @file WLFl.cpp
 * @brief Source of the ID for a sphere Worland(poly) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "QuICC/SpatialScheme/3D/WLFlBuilder.hpp"
#include "QuICC/Transform/SphereWorlandTransform.hpp"
#include "QuICC/Transform/SphereFftWorlandTransform.hpp"
#include "QuICC/Transform/Setup/Default.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/Equations/Tools/SHl.hpp"
#include "QuICC/Equations/Tools/SHlm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string WLFl::sTag = "WLFl";

   const std::string WLFl::sFormatted = "WLFl";

   const std::size_t WLFl::sId = Hasher::makeId(WLFl::sTag);

   WLFl::WLFl(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : xLF(formulation, purpose, WLFl::sId, WLFl::sTag, WLFl::sFormatted)
   {
      // Enable basic scheme features
      this->enable(Feature::FourierIndex2);
      this->enable(Feature::SpectralMatrix1D);
      this->enable(Feature::SpectralOrdering132);
   }

   void WLFl::setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type)
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

      // Set LF implementations
      xLF::setImplementation(type);
   }

   std::shared_ptr<IBuilder> WLFl::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<WLFlBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> WLFl::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            const auto& impl = this->mImplType.at(0);
            if(std::find(impl.begin(), impl.end(), Transform::Setup::Fft::id()) != impl.end())
            {
               auto spWT = std::make_shared<Transform::SphereFftWorlandTransform>();
               auto spST = std::dynamic_pointer_cast<Transform::SphereFftWorlandTransform::SetupType>(spSetup);
               if(!spST)
               {
                  throw std::logic_error("Incompatible transform setup for 1D (FFT) given");
               }
               spWT->init(spST);
               spTransform = spWT;
            } else if(std::find(impl.begin(), impl.end(), Transform::Setup::GaussianQuadrature::id()) != impl.end())
            {
               auto spWT = std::make_shared<Transform::SphereWorlandTransform>();
               auto spST = std::dynamic_pointer_cast<Transform::SphereWorlandTransform::SetupType>(spSetup);
               if(!spST)
               {
                  throw std::logic_error("Incompatible transform setup for 1D (Poly) (given");
               }
               spWT->init(spST);
               spTransform = spWT;
            }
            else
            {
               throw std::logic_error("Implementation type is not set properly for 1D");
            }
            break;
         }
         default:
            spTransform = xLF::createTransform(id, spSetup);
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
         default:
            spConv = xLF::createIndexConv(id);
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

} // SpatialScheme
} // QuICC
