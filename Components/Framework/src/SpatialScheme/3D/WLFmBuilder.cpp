/**
 * @file WLFmBuilder.cpp
 * @brief Source of the sphere Worland(poly) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral m ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFmBuilder.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Framework/MpiFramework.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/SpatialScheme/3D/WLFMesher.hpp"
#include "QuICC/Resolutions/Tools/RegularSHlIndexCounter.hpp"
#include "QuICC/Resolutions/Tools/RegularSHmIndexCounter.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   void WLFmBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
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

         const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA1D);

         // Initialise the ranks with local rank
         std::vector<std::set<int> >  ranks;
         ArrayI modes(tRes.dim<Dimensions::Data::DAT3D>());
         std::map<int, int>  mapModes;
         int k_ = 0;
         for(int k = 0; k < spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); ++k)
         {
            ranks.push_back(std::set<int>());
            if(k_ < tRes.dim<Dimensions::Data::DAT3D>() && k == tRes.idx<Dimensions::Data::DAT3D>(k_))
            {
               ranks.back().insert(QuICCEnv().id());
               modes(k_) = tRes.idx<Dimensions::Data::DAT3D>(k_);
               mapModes.insert(std::make_pair(tRes.idx<Dimensions::Data::DAT3D>(k_),k));
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

         MpiFramework::initSubComm(MpiFramework::SPECTRAL, tRes.dim<Dimensions::Data::DAT3D>());

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

            if(i_ < tRes.dim<Dimensions::Data::DAT3D>() && i == static_cast<size_t>(tRes.idx<Dimensions::Data::DAT3D>(i_)))
            {
               i_++;
            }
         }

         // Free communicator
         ierr = MPI_Comm_free(&comm);
         QuICCEnv().check(ierr, 825);
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
   }

   void WLFmBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      auto  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      auto  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      auto  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedTransformSetup WLFmBuilder::spSetup1D(SharedResolution spRes) const
   {
      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA1D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      Transform::SharedTransformSetup spSetup;

      const auto& opt = this->mOptions.at(0);

      // Poly algorithm setup
      if(std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end())
      {
         auto spPolySetup = std::make_shared<Transform::Poly::Setup>(size, specSize, this->purpose());
         spSetup = spPolySetup;
      }
      // FFT algorithm setup
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end())
      {
         auto spFftSetup = std::make_shared<Transform::Fft::Worland::Setup>(size, specSize, this->purpose());
         spSetup = spFftSetup;
      }

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         auto l = tRes.idx<Dimensions::Data::DAT3D>(i);
         auto nN = tRes.dim<Dimensions::Data::DATB1D>(0,i);

         spSetup->addIndex(l, tRes.dim<Dimensions::Data::DAT2D>(i), nN);
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::SharedTransformSetup WLFmBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(1);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA2D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::ALegendre::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(tRes.idx<Dimensions::Data::DAT3D>(i), tRes.dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::SharedTransformSetup WLFmBuilder::spSetup3D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(2);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA3D);

      // Get size of FFT transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += tRes.dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());
      spSetup->setBoxScale(1.0);
      spSetup->lock();

      return spSetup;
   }

   WLFmBuilder::WLFmBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegularSHmlBuilder(dim, purpose, options)
   {
   }

   bool WLFmBuilder::sameSpectralOrdering() const
   {
      return false;
   }

   void WLFmBuilder::setDimensions()
   {
      // Set default mesher
      auto m = std::make_shared<WLFMesher>(this->purpose());
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
