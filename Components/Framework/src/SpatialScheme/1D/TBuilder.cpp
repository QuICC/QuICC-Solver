/** 
 * @file TBuilder.cpp
 * @brief Source of the Chebyshev(FFT) scheme implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/1D/TBuilder.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void TBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);
   }

   Transform::Fft::Chebyshev::SharedSetup TBuilder::spSetup1D(SharedResolution spRes) const
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

   TBuilder::TBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular1DBuilder(dim, purpose, {})
   {
   }

   TBuilder::~TBuilder()
   {
   }

   void TBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(1);
      traSize(0) = this->mI + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::Fft::Tools::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::Fft::Tools::optimizeFft(nX);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nX, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);
   }

   void TBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);
   }

   void TBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);
   }

   void TBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);
   }

}
}
