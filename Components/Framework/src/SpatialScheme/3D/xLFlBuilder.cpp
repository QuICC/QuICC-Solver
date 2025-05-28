/** 
 * @file xLFlBuilder.cpp
 * @brief Source of the sphere Worland + Spherical Harmonics (Associated Legendre + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/xLFlBuilder.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace SpatialScheme {

   void xLFlBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for second transform
      auto  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      auto  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedTransformSetup xLFlBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(1);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA2D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::ALegendre::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(tRes.idx<Dimensions::Data::DAT3D>(i), tRes.dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::SharedTransformSetup xLFlBuilder::spSetup3D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(2);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA3D);

      // Get size of FFT transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += tRes.dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());
      spSetup->setBoxScale(1.0);
      spSetup->lock();

      return spSetup;
   }

   xLFlBuilder::xLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegularSHlBuilder(dim, purpose, options)
   {
   }

} // SpatialScheme
} // QuICC
