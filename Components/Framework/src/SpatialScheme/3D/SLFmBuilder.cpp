/** 
 * @file SLFmBuilder.cpp
 * @brief Source of the spherical Chebyshev(FFT) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with m spectral ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFmBuilder.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Resolutions/Tools/RegularSHlIndexCounter.hpp"
#include "QuICC/Resolutions/Tools/RegularSHmIndexCounter.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/SpatialScheme/3D/SLFMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   void SLFmBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      this->tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         // MPI error code
         int ierr;

         // Extract modes available on local rank
         std::vector<std::set<int> >  ranks;
         ArrayI modes(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
         std::map<int, int>  mapModes;
         int k_ = 0;
         for(int k = 0; k < spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); ++k)
         {
            ranks.push_back(std::set<int>());
            if(k_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && k == spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_))
            {
               ranks.back().insert(QuICCEnv().id());
               modes(k_) = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_);
               mapModes.insert(std::make_pair(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k_),k));
               k_++;
            }
         }

         // Share modes information with CPUs in first transform group
         int traId = 0;
         ArrayI tmp;
         for(int commCpu = 0; commCpu < MpiFramework::transformCpus(0).size(); ++commCpu)
         {
            ArrayI iData(2);
            iData.setConstant(-1);
            if(commCpu == MpiFramework::transformId(traId))
            {
               // Send the size and global rank
               iData(0) = modes.size();
               iData(1) = QuICCEnv().id();
               MpiFramework::syncTransform(traId);
               ierr = MPI_Bcast(iData.data(), iData.size(), MPI_INT, commCpu, MpiFramework::transformComm(traId));
               QuICCEnv().check(ierr, 815);

               // Send modes
               MpiFramework::syncTransform(traId);
               ierr = MPI_Bcast(modes.data(), modes.size(), MPI_INT, commCpu, MpiFramework::transformComm(traId));
               QuICCEnv().check(ierr, 817);
            } else
            {
               // Get size and global rank
               MpiFramework::syncTransform(traId);
               ierr = MPI_Bcast(iData.data(), iData.size(), MPI_INT, commCpu, MpiFramework::transformComm(traId));
               QuICCEnv().check(ierr, 818);

               // Receive modes
               tmp.resize(iData(0));
               MpiFramework::syncTransform(traId);
               ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, commCpu, MpiFramework::transformComm(traId));
               QuICCEnv().check(ierr, 820);

               // Expand modes to rank map
               std::map<int,int>::iterator mapIt;
               for(int i = 0; i < iData(0); i++)
               {
                  mapIt = mapModes.find(tmp(i));
                  if(mapIt != mapModes.end())
                  {
                     ranks.at(mapIt->second).insert(iData(1));
                  }
               }
            }
         }

         // Synchronize
         QuICCEnv().synchronize();

         // Initialize storage for SPECTRAL sub communicators
         MpiFramework::initSubComm(MpiFramework::SPECTRAL, spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

         std::set<int>  subRanks;
         int i_ = 0;
         for(size_t i = 0; i < ranks.size(); i++)
         {
            subRanks.clear();
            for(int cpu = 0; cpu < spRes->nCpu(); ++cpu)
            {
               int size;
               if(cpu == QuICCEnv().id())
               {
                  size = ranks.at(i).size();
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 821);

                  if(size > 0)
                  {
                     tmp.resize(size);
                     int j = 0;
                     for(auto sIt = ranks.at(i).begin(); sIt != ranks.at(i).end(); ++sIt)
                     {
                        tmp(j) = *sIt;
                        ++j;
                        subRanks.insert(*sIt);
                     }
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(tmp.data(), size, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 822);
                  }
               } else
               {
                  // Get size
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 823);

                  // Receive ranks
                  if(size > 0)
                  {
                     tmp.resize(size);
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 824);

                     for(int j = 0; j < size; ++j)
                     {
                        subRanks.insert(tmp(j));
                     }
                  }
               }
            }

            MpiFramework::setSubComm(MpiFramework::SPECTRAL, i_, subRanks);

            if(i_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && i == static_cast<size_t>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i_)))
            {
               i_++;
            }
         }

         // Synchronize
         QuICCEnv().synchronize();
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void SLFmBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Chebyshev::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Poly::ALegendre::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Chebyshev::SharedSetup SLFmBuilder::spSetup1D(SharedResolution spRes) const
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

   Transform::Poly::ALegendre::SharedSetup SLFmBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::ALegendre::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT3D>(i), spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup SLFmBuilder::spSetup3D(SharedResolution spRes) const
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
      spSetup->setBoxScale(1.0);
      spSetup->lock();

      return spSetup;
   }

   SLFmBuilder::SLFmBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegularSHmlBuilder(dim, purpose, options)
   {
      this->mOptions.at(0).push_back(Transform::Setup::Uniform::id());
   }

   bool SLFmBuilder::sameSpectralOrdering() const
   {
      return false;
   }

   void SLFmBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<SLFMesher>(this->purpose());
      this->setMesher(m, false);
      // ... initialize mesher
      std::vector<int> d = {this->mI, this->mL, this->mM};
      this->mesher().init(d, this->mOptions);

      // Set dimensions using mesher
      I3DBuilder::setDimensions();

      //
      // Change order for 2D, 3D in spectral space
      //

      // Initialise second dimension of first transform
      this->setDimension(this->mesher().nSpec2D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(this->mesher().nSpec3D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT3D);
   }

} // SpatialScheme
} // QuICC
