/** 
 * @file WFTBuilder.cpp
 * @brief Source of the cylindrical Worland(poly) + Fourier + Chebyshev(FFT) scheme implementation
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/WFTBuilder.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Framework/MpiFramework.hpp"

namespace QuICC {

namespace SpatialScheme {

   void WFTBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      this->tuneMpiResolution(descr);
      
      // Create spectral space sub communicators
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         // MPI error code
         int ierr;

         // Get world group
         MPI_Group world;
         MPI_Group group;
         ierr = MPI_Comm_group(MPI_COMM_WORLD, &world);
         QuICCEnv().check(ierr, 811);

         // Create minimial MPI group
         ierr = MPI_Group_incl(world, MpiFramework::transformCpus(0).size(), MpiFramework::transformCpus(0).data(), &group);
         QuICCEnv().check(ierr, 812);

         // Create minimial MPI communicator
         MPI_Comm comm;
         ierr = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
         QuICCEnv().check(ierr, 813);

         // Initialise the ranks with local rank
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

         // Loop over all cpus
         int commId;
         int globalCpu = QuICCEnv().id();
         ierr = MPI_Comm_rank(comm, &commId); 
         QuICCEnv().check(ierr, 814);
         ArrayI tmp;
         for(int commCpu = 0; commCpu < MpiFramework::transformCpus(0).size(); ++commCpu)
         {
            int size;
            if(commCpu == commId)
            {
               // Send the size
               size = modes.size();
               ierr = MPI_Bcast(&size, 1, MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 815);
               MPI_Barrier(comm);

               // Send global CPU rank 
               globalCpu = QuICCEnv().id();
               ierr = MPI_Bcast(&globalCpu, 1, MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 816);
               MPI_Barrier(comm);

               // Send modes
               ierr = MPI_Bcast(modes.data(), modes.size(), MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 817);
               MPI_Barrier(comm);
            } else
            {
               // Get size
               ierr = MPI_Bcast(&size, 1, MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 818);
               MPI_Barrier(comm);

               // Get global CPU rank 
               ierr = MPI_Bcast(&globalCpu, 1, MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 819);
               MPI_Barrier(comm);

               // Receive modes
               tmp.resize(size);
               ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, commCpu, comm);
               QuICCEnv().check(ierr, 820);
               MPI_Barrier(comm);

               std::map<int,int>::iterator mapIt;
               for(int i = 0; i < size; i++)
               {
                  mapIt = mapModes.find(tmp(i));
                  if(mapIt != mapModes.end())
                  {
                     ranks.at(mapIt->second).insert(globalCpu);
                  }
               }
            }

            // Synchronize
            QuICCEnv().synchronize();
         }

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
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 821);
                  QuICCEnv().synchronize();

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
                     ierr = MPI_Bcast(tmp.data(), size, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 822);
                     QuICCEnv().synchronize();
                  }
               } else
               {
                  // Get size
                  ierr = MPI_Bcast(&size, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 823);
                  QuICCEnv().synchronize();

                  // Receive ranks
                  if(size > 0)
                  {
                     tmp.resize(size);
                     ierr = MPI_Bcast(tmp.data(), tmp.size(), MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 824);
                     QuICCEnv().synchronize();

                     for(int j = 0; j < size; ++j)
                     {
                        subRanks.insert(tmp(j));
                     }
                  }
               }

               // Synchronize
               QuICCEnv().synchronize();
            }

            MpiFramework::setSubComm(MpiFramework::SPECTRAL, i_, subRanks);

            if(i_ < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>() && i == static_cast<size_t>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i_)))
            {
               i_++;
            }
         }
         
         // Free communicator
         ierr = MPI_Comm_free(&comm);
         QuICCEnv().check(ierr, 825);
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void WFTBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Poly::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Chebyshev::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Poly::SharedSetup WFTBuilder::spSetup1D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT3D>(i), spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup WFTBuilder::spSetup2D(SharedResolution spRes) const
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

   Transform::Fft::Chebyshev::SharedSetup WFTBuilder::spSetup3D(SharedResolution spRes) const
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

   WFTBuilder::WFTBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegular3DBuilder(dim, purpose)
   {
   }

   WFTBuilder::~WFTBuilder()
   {
   }

   void WFTBuilder::setDimensions()
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

      // Get dealiased Worland size
      int nR = Transform::Poly::Tools::dealias(this->mI+this->mJ/2+8);

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
      this->setDimension(this->mI+1, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

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

   void WFTBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void WFTBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void WFTBuilder::setMemoryScore()
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
