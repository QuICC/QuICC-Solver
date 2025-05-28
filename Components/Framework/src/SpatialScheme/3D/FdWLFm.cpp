/** 
 * @file FdWLFm.cpp
 * @brief Source of ID of the Finite Differences sphere + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral m ordering
 */

// System includes
//

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/3D/FdWLFm.hpp"
#include "QuICC/SpatialScheme/3D/FdWLFmBuilder.hpp"
#include "QuICC/Transform/SphereFiniteDiffTransform.hpp"
#include "QuICC/Transform/Setup/Default.hpp"
#include "QuICC/Transform/Setup/FiniteDiff.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/Communicators/Converters/SHm2lIndexConv.hpp"
#include "QuICC/Equations/Tools/SHm.hpp"

namespace QuICC {

namespace SpatialScheme {

   const std::string FdWLFm::sTag = "FdWLFm";

   const std::string FdWLFm::sFormatted = "FdWLFm";

   const std::size_t FdWLFm::sId = Hasher::makeId(FdWLFm::sTag);

   FdWLFm::FdWLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose)
      : xLF(formulation, purpose, FdWLFm::sId, FdWLFm::sTag, FdWLFm::sFormatted)
   {
      // Enable basic scheme features
      this->enable(Feature::FourierIndex3);
      this->enable(Feature::SpectralMatrix2D);
      this->enable(Feature::SpectralOrdering123);
   }

   void FdWLFm::setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type)
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
         mOpt1D.push_back(Transform::Setup::FiniteDiff::id());
         mOpt1D.push_back(Transform::Setup::Uniform::id());
      }
      else
      {
         mOpt1D = opt1D;
      }


      // Set LF implementations
      xLF::setImplementation(type);
   }

   std::shared_ptr<IBuilder> FdWLFm::createBuilder(ArrayI& dim, const bool needInterpretation) const
   {
      auto spBuilder = makeBuilder<FdWLFmBuilder>(dim, this->purpose(), needInterpretation, this->mImplType, this->mspCustomMesher);

      return spBuilder;
   }

   std::shared_ptr<Transform::ITransform> FdWLFm::createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const
   {
      std::shared_ptr<Transform::ITransform> spTransform;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
         {
            const auto& impl = this->mImplType.at(0);
            if(std::find(impl.begin(), impl.end(), Transform::Setup::FiniteDiff::id()) != impl.end())
            {
               auto spWT = std::make_shared<Transform::SphereFiniteDiffTransform>();
               auto spST = std::dynamic_pointer_cast<Transform::SphereFiniteDiffTransform::SetupType>(spSetup);
               if(!spST)
               {
                  throw std::logic_error("Incompatible transform setup for 1D FiniteDiff given");
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

   std::shared_ptr<Parallel::IIndexConv> FdWLFm::createIndexConv(const Dimensions::Transform::Id id) const
   {
      std::shared_ptr<Parallel::IIndexConv> spConv;

      switch(id)
      {
         case Dimensions::Transform::TRA1D:
            spConv = std::make_shared<Parallel::SHm2lIndexConv>();
            break;
         default:
            spConv = xLF::createIndexConv(id);
      }

      return spConv;
   }

   std::shared_ptr<Equations::Tools::ICoupling> FdWLFm::createCouplingTools(const Equations::CouplingIndexType indexType) const
   {
      if(indexType != Equations::CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         throw std::logic_error("Unknown coupling tools");
      }

      auto spCoupling = std::make_shared<Equations::Tools::SHm>();
      return spCoupling;
   }

} // SpatialScheme
} // QuICC
