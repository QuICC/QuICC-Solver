/** 
 * @file TTBuilder.cpp
 * @brief Source of the Chebyshev(FFT) + Chebyshev(FFT) scheme implementation
 */

// System includes
//
#include <set>
#include <vector>

// Project includes
//
#include "QuICC/SpatialScheme/2D/TTBuilder.hpp"
#include "QuICC/SpatialScheme/2D/TTMesher.hpp"
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

   TTBuilder::TTBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegular2DBuilder(dim, purpose, options)
   {
   }

   void TTBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<TTMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mJ};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      IRegular2DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
