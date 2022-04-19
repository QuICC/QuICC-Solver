/** 
 * @file Resolution.cpp
 * @brief Source of the resolution object for several CPUs
 */

// Configuration includes
//

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Resolutions/Resolution.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

   Resolution::Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim, const ArrayI& transDim)
      : mLocalId(-1), mCores(coreRes)
   {
      this->setLocalId();

      // Assert simulation dimensions are the same as cpu dimensions
      assert(simDim.size() == this->cpu()->nDim());

      // Init the simulation resolution
      this->initSimResolution(simDim, transDim);
   }

   Resolution::~Resolution()
   {
   }

   void Resolution::setLocalId()
   {
      // Loop over core resolutions
      for(std::size_t i = 0; i < this->mCores.size(); i++)
      {
         if(!this->mCores.at(i)->isCleared())
         {
            if(this->mLocalId != -1)
            {
               throw std::logic_error("Multiple local ID obtained from core resolutions");
            }
            this->mLocalId = i;
         }
      }

      if(this->mLocalId == -1)
      {
         throw std::logic_error("Local ID could not be obtained from core resolutions");
      }
   }

   void Resolution::initSimResolution(const ArrayI& simDim, const ArrayI& transDim)
   {
      // Get number of dimensions
      int nDim = this->cpu()->nDim();

      // Create storage for the physical dimensions
      ArrayI phys(nDim);
      for(int i = 0; i < nDim; i++)
      {
         phys(i) = this->cpu()->dim(static_cast<Dimensions::Transform::Id>(i))->dim<Dimensions::Data::DATF1D>();
      }

      ArrayI spec = simDim;
      spec.array() += 1;

      // Create shared pointer of simulation resolution
      this->mspSim = std::make_shared<SimulationResolution>(phys, spec, transDim);
   }

   int Resolution::nCpu() const
   {
      return this->mCores.size();
   }

   void Resolution::setBoxScale(const Array& boxScale)
   {
      this->mspSim->setBoxScale(boxScale);

      // Set the boxscale on the 
      for(auto it = this->mTSetups.begin(); it != this->mTSetups.end(); ++it)
      {
         it->second->setBoxScale(this->sim().boxScale(static_cast<Dimensions::Simulation::Id>(static_cast<int>(it->first))));
      }
      
   }

   void Resolution::setSpatialScheme(SpatialScheme::SharedISpatialScheme spScheme)
   {
      this->mspSim->setSpatialScheme(spScheme);
   }

   void Resolution::setIndexCounter(SharedIndexCounter spCounter)
   {
      this->mspCounter = spCounter;
   }

   const IndexCounter& Resolution::counter() const
   {
      // Safety assert
      assert(this->mspCounter);

      return *this->mspCounter;
   }

   SharedCSimulationResolution Resolution::spSim() const
   {
      // Safety assert
      assert(this->mspSim);

      return this->mspSim;
   }

   const SimulationResolution& Resolution::sim() const
   {
      // Safety assert
      assert(this->mspSim);

      return *this->mspSim;
   }

   SharedCCoreResolution Resolution::spCpu() const
   {
      // Safety assert
      assert(this->mLocalId >= 0);
      assert(this->mCores.size() > static_cast<size_t>(this->mLocalId));

      return this->mCores.at(this->mLocalId);
   }

   SharedCCoreResolution Resolution::cpu() const
   {
      // Safety assert
      assert(this->mLocalId >= 0);
      assert(this->mCores.size() > static_cast<size_t>(this->mLocalId));

      return this->mCores.at(this->mLocalId);
   }

   SharedCCoreResolution Resolution::cpu(const int id) const
   {
      // Check sizes
      assert(id >= 0);
      assert(static_cast<size_t>(id) < this->mCores.size());
      
      return this->mCores.at(id);
   }

   void Resolution::addTransformSetup(const Dimensions::Transform::Id id, Transform::SharedTransformSetup spSetup)
   {
      this->mTSetups.insert(std::make_pair(id, spSetup));
   }

   Transform::SharedTransformSetup  Resolution::spTransformSetup(const Dimensions::Transform::Id id) const
   {
      // Assert for correct size
      assert(this->mTSetups.size() > static_cast<size_t>(id));

      return this->mTSetups.find(id)->second;
   }

   std::shared_ptr<Datatypes::ScalarFieldSetup> Resolution::spFwdSetup(const Dimensions::Transform::Id id) const
   {
      // Get forward dimensions
      auto spDim1D = std::make_shared<ArrayI>(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < spDim1D->size(); ++i)
      {
         (*spDim1D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DATF1D>(i);
      }

      // Get 2D dimensions
      auto spDim2D = std::make_shared<ArrayI>(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DAT2D>(i);
      }

      return std::make_shared<Datatypes::ScalarFieldSetup>(spDim1D, spDim2D, this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
   }

   std::shared_ptr<Datatypes::ScalarFieldSetup> Resolution::spBwdSetup(const Dimensions::Transform::Id id) const
   {
      // Get backward dimensions
      auto spDim1D = std::make_shared<ArrayI>(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < spDim1D->size(); ++i)
      {
         (*spDim1D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DATB1D>(i);
      }

      // Get 2D dimensions
      auto spDim2D = std::make_shared<ArrayI>(this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(id)->dim<Dimensions::Data::DAT2D>(i);
      }

      return std::make_shared<Datatypes::ScalarFieldSetup>(spDim1D, spDim2D, this->cpu()->dim(id)->dim<Dimensions::Data::DAT3D>());
   }

   std::shared_ptr<Datatypes::ScalarFieldSetup> Resolution::spSpectralSetup() const
   {
      // Get backward dimensions
      auto spDim1D = std::make_shared<ArrayI>(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
      spDim1D->setConstant(this->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM));
      for(int i = 0; i < spDim1D->size(); ++i)
      {
         (*spDim1D)(i) = this->counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i));
      }

      // Get 2D dimensions
      auto spDim2D = std::make_shared<ArrayI>(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < spDim2D->size(); ++i)
      {
         (*spDim2D)(i) = this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return std::make_shared<Datatypes::ScalarFieldSetup>(spDim1D, spDim2D, this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
   }

   std::shared_ptr<Datatypes::ScalarFieldSetup> Resolution::spPhysicalSetup() const
   {
      return this->spFwdSetup(static_cast<Dimensions::Transform::Id>(this->cpu()->nDim()-1));
   }

   void Resolution::buildRestriction(std::vector<int>& rSlow, std::vector<std::vector<int> >& rMiddle, const int k)
   {
      // Make sure vectors are empty
      rSlow.clear();
      rMiddle.clear();

      #ifdef QUICC_MPISPSOLVE
         // Three dimensional matrices
         if(this->sim().ss().has(SpatialScheme::Feature::SpectralMatrix3D))
         {
            int k_;
            int j_;
            rSlow.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
            rMiddle.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());

            for(int k=0; k < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Extract "physical" index of slow data dimension
               k_ = this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

               rSlow.push_back(k_);
               rMiddle.push_back(std::vector<int>());
               rMiddle.back().reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k));

               // Loop over middle data dimension
               for(int j=0; j < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  // Extract "physical" index of middle data dimension
                  j_ = this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);

                  rMiddle.back().push_back(j_);
               }
            }

         // Two dimensional matrices
         } else if(this->sim().ss().dimension() == 3 && this->sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
         {
            // Make sure middle is empty
            rMiddle.clear();

            rSlow.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k));

            for(int j = 0; j < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               rSlow.push_back(this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k));
            }

         // Two dimensional matrices in 2D simulation
         } else if(this->sim().ss().dimension() == 2 && this->sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
         {
            // Make sure middle is empty
            rMiddle.clear();

            rSlow.reserve(this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k));

            for(int j = 0; j < this->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               rSlow.push_back(this->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k));
            }

         } else
         {
            throw std::logic_error("Tried to setup MPI sparse solver for unimplemented scheme.");
         }
      #else
         throw std::logic_error("buildRestriction was called without MPI sparse solver");
      #endif //QUICC_MPISPSOLVE
   }

}
