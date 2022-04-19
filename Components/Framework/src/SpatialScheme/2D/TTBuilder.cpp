/** 
 * @file TTBuilder.cpp
 * @brief Source of the Chebyshev(FFT) + Chebyshev(FFT) scheme implementation
 */

// System includes
//

// External includes
//
#include <set>
#include <vector>

// Class include
//
#include "QuICC/SpatialScheme/2D/TTBuilder.hpp"

// Project includes
//
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Resolutions/Tools/RegularIndexCounter.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void TTBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      this->tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         // Initialise the ranks with local rank
         std::set<int>  ranks;
         for(int cpu = 0; cpu < spRes->nCpu(); ++cpu)
         {
            ranks.insert(cpu);
         }

         MpiFramework::initSubComm(MpiFramework::SPECTRAL, 1);

         MpiFramework::setSubComm(MpiFramework::SPECTRAL, 0, ranks);
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void TTBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Chebyshev::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);
   }

   Transform::Fft::Chebyshev::SharedSetup TTBuilder::spSetup1D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>();

      auto spSetup = std::make_shared<Transform::Fft::Chebyshev::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Chebyshev::SharedSetup TTBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>();

      auto spSetup = std::make_shared<Transform::Fft::Chebyshev::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   TTBuilder::TTBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular2DBuilder(dim, purpose)
   {
   }

   TTBuilder::~TTBuilder()
   {
   }

   void TTBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(2);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mJ + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get standard dealiased FFT size
      int nX = Transform::Fft::Tools::dealiasCosFft(this->mI+1);
      // Check for optimised FFT sizes
      nX = Transform::Fft::Tools::optimizeFft(nX);

      // Get standard dealiased FFT size
      int nY = Transform::Fft::Tools::dealiasCosFft(this->mJ+1);
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
      this->setDimension(nY, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nX, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);
   }

   void TTBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);
   }

   void TTBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);
   }

   void TTBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);
   }

   bool TTBuilder::applicable() const
   {
      return true;
   }

}
}
