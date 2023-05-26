/** 
 * @file TFFBuilder.cpp
 * @brief Source of the Chebyshev(FFT) + Fourier + Fourier scheme implementation
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/TFFBuilder.hpp"
#include "QuICC/Resolutions/Tools/DoublePeriodicIndexCounter.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   void TFFBuilder::interpretConfigDimensions(ArrayI& rDim)
   {
      rDim(1) = 2*rDim(1);
   }

   void TFFBuilder::addIndexCounter(SharedResolution spRes)
   {
      auto spCounter = std::make_shared<DoublePeriodicIndexCounter>(spRes->spSim(), spRes->spCpu());

      spRes->setIndexCounter(spCounter);
   }

   void TFFBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Fourier::Complex::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Chebyshev::SharedSetup TFFBuilder::spSetup1D(SharedResolution spRes) const
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

      // Check for mean
      MatrixI idBlocks;
      if(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0 && spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0) == 0)
      {
         spSetup->addIndex(0,1);
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Complex::SharedSetup TFFBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Complex::Setup>(size, specSize, this->purpose());

      // Get number of transforms
      std::vector<std::pair<int,int> > idPairs;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         // Store 3D ID and block size
         spSetup->addIndex(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT3D>(i), spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup TFFBuilder::spSetup3D(SharedResolution spRes) const
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

   TFFBuilder::TFFBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular3DBuilder(dim, purpose, {})
   {
      this->mOptions.at(0).push_back(Transform::Setup::Uniform::id());
   }

   void TFFBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      /// \mhdBug possible aliasing issue, check for correct dealiased FFT
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mJ + 1;
      traSize(2) = this->mK + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased cosine FFT size
      int nX = Transform::Fft::Tools::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::Fft::Tools::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::Fft::Tools::dealiasFft(this->mJ+1);
      // Check for optimised FFT sizes
      nY = Transform::Fft::Tools::optimizeFft(nY);

      // Get standard dealiased FFT size
      int nZ = Transform::Fft::Tools::dealiasMixedFft(this->mK+1);
      // Check for optimised FFT sizes
      nZ = Transform::Fft::Tools::optimizeFft(nZ);

      //
      // Initialise spectral space
      //

      // Initialise forward dimension of first transform
      this->setDimension(traSize(0), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(traSize(0), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT3D);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nZ, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nZ/2 + 1, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nY, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nX, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void TFFBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void TFFBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void TFFBuilder::setMemoryScore()
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
