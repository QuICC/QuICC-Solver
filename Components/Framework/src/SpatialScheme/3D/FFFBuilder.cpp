/** 
 * @file FFFBuilder.cpp
 * @brief Source of the Fourier + Fourier + Fourier scheme implementation
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/FFFBuilder.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/SpatialScheme/3D/FFFMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   void FFFBuilder::interpretConfigDimensions(ArrayI& rDim)
   {
      rDim(0) = 2*rDim(0);
      rDim(1) = 2*rDim(1);
   }

   void FFFBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Fourier::Complex::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Fourier::Complex::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Fourier::Complex::SharedSetup FFFBuilder::spSetup1D(SharedResolution spRes) const
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

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Complex::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Complex::SharedSetup FFFBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Complex::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup FFFBuilder::spSetup3D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   FFFBuilder::FFFBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegular3DBuilder(dim, purpose, options)
   {
   }

   void FFFBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<FFFMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mJ, this->mK};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      I3DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
