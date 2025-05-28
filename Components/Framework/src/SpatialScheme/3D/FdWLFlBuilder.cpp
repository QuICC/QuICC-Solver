/** 
 * @file FdWLFlBuilder.cpp
 * @brief Source of the Finite Differnces sphere + Spherical Harmonics (Associated Legendre + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/FdWLFlBuilder.hpp"
#include "QuICC/Transform/FiniteDiff/Setup.hpp"
#include "QuICC/Transform/Setup/FiniteDiff.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SpatialScheme/3D/FdWLFMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   void FdWLFlBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      auto  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Setup 2D/3D
      xLFlBuilder::addTransformSetups(spRes);
   }

   Transform::SharedTransformSetup FdWLFlBuilder::spSetup1D(SharedResolution spRes) const
   {
      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA1D);

      // Get size of the transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      Transform::SharedTransformSetup spSetup;

      const auto& opt = this->mOptions.at(0);

      // Finite Differences algorithm setup
      if(std::find(opt.begin(), opt.end(), Transform::Setup::FiniteDiff::id()) != opt.end())
      {
         auto spFdSetup = std::make_shared<Transform::FiniteDiff::Setup>(specSize, this->purpose());
         spSetup = spFdSetup;
      }
      else
      {
         throw std::logic_error("Unknown finite differences algorithm");
      }

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         auto l = tRes.idx<Dimensions::Data::DAT3D>(i);
         auto nN = tRes.dim<Dimensions::Data::DATB1D>(0,i);

         spSetup->addIndex(l, tRes.dim<Dimensions::Data::DAT2D>(i), nN);
      }

      spSetup->lock();

      return spSetup;
   }

   FdWLFlBuilder::FdWLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : xLFlBuilder(dim, purpose, options)
   {
   }

   void FdWLFlBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<FdWLFMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mL, this->mM};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      I3DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
