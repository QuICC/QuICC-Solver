/** 
 * @file FFBuilder.cpp
 * @brief Source of the Fourier + Fourier scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/2D/FFBuilder.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void FFBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Fourier::Complex::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);
   }

   Transform::Fft::Fourier::Complex::SharedSetup FFBuilder::spSetup1D(SharedResolution spRes) const
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

   Transform::Fft::Fourier::Mixed::SharedSetup FFBuilder::spSetup2D(SharedResolution spRes) const
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

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   FFBuilder::FFBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular2DBuilder(dim, purpose, {})
   {
   }

   FFBuilder::~FFBuilder()
   {
   }

   void FFBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(2);
      traSize(0) = this->mI + 1;
      traSize(1) = 2*(this->mJ + 1);
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::Fft::Tools::dealiasFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::Fft::Tools::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::Fft::Tools::dealiasMixedFft(this->mJ+1);
      // Check for optimised FFT sizes
      nY = Transform::Fft::Tools::optimizeFft(nY);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nY/2 + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);
   }

   void FFBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);
   }

   void FFBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);
   }

   void FFBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);
   }


}
}
