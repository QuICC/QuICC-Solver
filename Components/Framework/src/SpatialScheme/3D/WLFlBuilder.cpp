/** 
 * @file WLFlBuilder.cpp
 * @brief Source of the sphere Worland + Spherical Harmonics (Associated Legendre + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFlBuilder.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SpatialScheme/3D/WLFMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   void WLFlBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      auto  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Setup 2D/3D
      xLFlBuilder::addTransformSetups(spRes);
   }

   Transform::SharedTransformSetup WLFlBuilder::spSetup1D(SharedResolution spRes) const
   {
      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA1D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      Transform::SharedTransformSetup spSetup;

      const auto& opt = this->mOptions.at(0);

      // Poly algorithm setup
      if(std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end())
      {
         auto spPolySetup = std::make_shared<Transform::Poly::Setup>(size, specSize, this->purpose());
         spSetup = spPolySetup;
      }
      // FFT algorithm setup
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end())
      {
         auto spFftSetup = std::make_shared<Transform::Fft::Worland::Setup>(size, specSize, this->purpose());
         spSetup = spFftSetup;
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

   WLFlBuilder::WLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : xLFlBuilder(dim, purpose, options)
   {
   }

   void WLFlBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<WLFMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mL, this->mM};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      I3DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
