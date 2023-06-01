/** 
 * @file TTTBuilder.cpp
 * @brief Source of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme implementation
 */

// System includes
//
#include <set>
#include <vector>

// Project includes
//
#include "QuICC/SpatialScheme/3D/TTTBuilder.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Tools/RegularIndexCounter.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/SpatialScheme/3D/TTTMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   void TTTBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      // Create single rank communicator
      #ifdef QUICC_MPI
         #pragma message("TTTBuilder seems to be broken in MPI mode!")
         std::vector<int>  ranks;

         // Make single core communicator unless MPI sparse solver is used
         #ifndef QUICC_MPISPSOLVE
            ranks.push_back(QuICCEnv().id());
         #endif //QUICC_MPISPSOLVE

         //MpiFramework::setSpectralComm(ranks);
      #endif //QUICC_MPI
   }

   void TTTBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Chebyshev::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Chebyshev::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Chebyshev::SharedSetup TTTBuilder::spSetup1D(SharedResolution spRes) const
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

   Transform::Fft::Chebyshev::SharedSetup TTTBuilder::spSetup2D(SharedResolution spRes) const
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

      auto spSetup = std::make_shared<Transform::Fft::Chebyshev::Setup>(size, blockSize, specSize, this->purpose());

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Chebyshev::SharedSetup TTTBuilder::spSetup3D(SharedResolution spRes) const
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

   TTTBuilder::TTTBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegular3DBuilder(dim, purpose, options)
   {
   }

   void TTTBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<TTTMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mJ, this->mK};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      I3DBuilder::setDimensions();
   }

} // SpatialScheme
} // QuICC
