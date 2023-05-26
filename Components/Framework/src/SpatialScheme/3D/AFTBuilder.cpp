/** 
 * @file AFTBuilder.cpp
 * @brief Source of the annulus Chebyshev(FFT) + Fourier + Chebyshev(FFT) scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/AFTBuilder.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void AFTBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Chebyshev::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Chebyshev::SharedSetup AFTBuilder::spSetup1D(SharedResolution spRes) const
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

      auto spSetup = std::make_shared<Transform::Fft::Chebyshev::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup AFTBuilder::spSetup2D(SharedResolution spRes) const
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

   Transform::Fft::Chebyshev::SharedSetup AFTBuilder::spSetup3D(SharedResolution spRes) const
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

      auto spSetup = std::make_shared<Transform::Fft::Chebyshev::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   AFTBuilder::AFTBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular3DBuilder(dim, purpose, {})
   {
   }

   AFTBuilder::~AFTBuilder()
   {
   }

   void AFTBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mJ + 1;
      traSize(2) = this->mK + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nR = Transform::Fft::Tools::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nR = Transform::Fft::Tools::optimizeFft(nR);

      // Get mixed dealiased FFT size
      int nTh = Transform::Fft::Tools::dealiasMixedFft(this->mJ+1);
      // Check for optimised FFT sizes
      nTh = Transform::Fft::Tools::optimizeFft(nTh);

      // Get standard dealiased FFT size
      int nZ = Transform::Fft::Tools::dealiasCosFft(this->mK+1);
      // Check for optimised FFT sizes
      nZ = Transform::Fft::Tools::optimizeFft(nZ);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nTh, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nTh/2 + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nR, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nTh, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nR, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void AFTBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void AFTBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void AFTBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

}
}
