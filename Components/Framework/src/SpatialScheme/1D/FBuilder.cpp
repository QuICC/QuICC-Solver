/** 
 * @file FBuilder.cpp
 * @brief Source of the Fourier scheme implementation
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/1D/FBuilder.hpp"
#include "QuICC/SpatialScheme/1D/FMesher.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void FBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);
   }

   Transform::Fft::Fourier::Mixed::SharedSetup FBuilder::spSetup1D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();
      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   FBuilder::FBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegular1DBuilder(dim, purpose, options)
   {
   }

   void FBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<FMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      IRegular1DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
