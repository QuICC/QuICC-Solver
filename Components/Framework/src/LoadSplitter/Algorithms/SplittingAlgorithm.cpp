/**
 * @file SplittingAlgorithm.cpp
 * @brief Source of the base of the implementation of the load splitting algorithms
 */

// System includes
//
#include <algorithm>
#include <map>
#include <stdexcept>

#ifdef QUICC_MPI
#include <mpi.h>
#endif

// External includes
//

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithmDetails.hpp"
#include "Profiler/Interface.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Parallel {

   SplittingAlgorithm::SplittingAlgorithm(const int id, const int nCpu, const ArrayI& dim, const Splitting::Algorithms::Id algo)
      : mAlgo(algo), mGrouper(Splitting::Groupers::EQUATION), mId(id), mNCpu(nCpu), mDims(dim.size()), mSimDim(dim)
   {
   }

   int SplittingAlgorithm::id() const
   {
      return this->mId;
   }

   int SplittingAlgorithm::nCpu() const
   {
      return this->mNCpu;
   }

   int SplittingAlgorithm::dims() const
   {
      return this->mDims;
   }

   int SplittingAlgorithm::factor(const int i) const
   {
      // Assert on index of requested factor
      assert(i < this->mFactors.size());

      return this->mFactors(i);
   }

   void SplittingAlgorithm::useFactorization(const std::list<int>& f)
   {
      // Use imposed factors is not empty
      if(f.size() > 0)
      {
         // Check size
         if(f.size() % this->mFactors.size() != 0)
         {
            throw std::logic_error("List of factors is not compatible");
         }

         this->mNCpuFactors.clear();
         this->mNCpuFactors = f;

         // Check factorization
         SplittingTools::filterFactors(this->mNCpuFactors, this->mFactors.size(), this->nCpu(), false);
      }
   }

   const ArrayI& SplittingAlgorithm::factors() const
   {
      return this->mFactors;
   }

   int SplittingAlgorithm::maxFactor() const
   {
      return this->mFactors.maxCoeff();
   }

   void SplittingAlgorithm::setScheme(SpatialScheme::SharedIBuilder spBuilder)
   {
      this->mspScheme = spBuilder;
   }

   std::pair<int, std::pair<SharedResolution, SplittingDescription> > SplittingAlgorithm::scoreSplitting(const Splitting::Groupers::Id grp)
   {
      StageTimer stage;
      std::stringstream ss;
      ss << this->factor(0);
      for(int i = 1; i < this->mFactors.size(); i++)
      {
         ss << " x " << this->factor(i);
      }
      stage.start("computing load splitting " + ss.str(), 1);

      // Storage for all the shared core resolutions
      std::vector<SharedCoreResolution>  coreRes;

      // Storage for all the shared transform resolutions
      std::vector<SharedTransformResolution>  transformRes;

      // Initialise description
      SplittingDescription descr;

      // Load splitting might fail
      int status = 0;

      // Loop over all CPUs
      for(int id = 0; id < this->nCpu(); id++)
      {
         // Clear content of the transform resolutions
         transformRes.clear();

         // Loop over all dimensions
         for(int j = 0; j < this->dims(); j++)
         {
            SharedTransformResolution  spTRes = this->splitDimension(static_cast<Dimensions::Transform::Id>(j), id, status);

            QuICCEnv().synchronize();
            #ifdef QUICC_MPI
               MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            #endif //QUICC_MPI

            // Splitting fail, abort
            if(status != 0)
            {
               break;
            }

            // Add stage to description
            descr.addStage(j, spTRes, id);

            // Clear unused indexes for remote resolutions
            if(id != this->id())
            {
               spTRes->clearIndexes();
            }

            transformRes.push_back(spTRes);
         }

         // Splitting fail, abort
         if(status != 0)
         {
            break;
         }

         // Create spectral resolution
         SharedTransformResolution  spSpectralRes = this->splitDimension(Dimensions::Transform::SPECTRAL, id, status);

         QuICCEnv().synchronize();
         #ifdef QUICC_MPI
            MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         #endif //QUICC_MPI

         // Splitting fail, abort
         if(status != 0)
         {
            break;
         }

         // Add stage to description
         descr.addStage(static_cast<int>(Dimensions::Transform::SPECTRAL), spSpectralRes, id);

         // Clear unused indexes for remote resolutions
         if(id != this->id())
         {
            spSpectralRes->clearIndexes();
         }

         // Create new shared core resolution
         coreRes.push_back(std::make_shared<CoreResolution>(transformRes, spSpectralRes));
      }

      stage.done();
      stage.start("creating resolution for " + ss.str(), 1);

      SharedResolution  spRes;

      // Splitting was successful
      Array score = Array::Constant(4,1.0);
      if(status == 0)
      {
         // Create shared resolution
         ArrayI transDim = this->mspScheme->getTransformSpace();
         spRes = std::make_shared<Resolution>(coreRes, this->mSimDim, transDim);

         // Add index counter to resolution
         this->mspScheme->addIndexCounter(spRes);

         // Add the transform setups to the resolution
         this->mspScheme->addTransformSetups(spRes);

         // Compute the score of the obtained resolution
         score = this->computeScore(spRes, grp);
      } else
      {
         // Set large negative score (splitting is unusable)
         score(0) = -9999;
      }

      // Create splitting description
      descr.algorithm = this->mAlgo;
      descr.grouper = this->mGrouper;
      descr.dims = this->mDims;
      descr.factors = this->mFactors;
      descr.score = score;
      descr.structure = std::map<Dimensions::Transform::Id, std::multimap<int,int> >();

      stage.done();

      // Return combination of score and shared resolution/description
      return std::make_pair(static_cast<int>(score.prod()), std::make_pair(spRes,descr));
   }

   void SplittingAlgorithm::initFactors(const int nFactors)
   {
      // Initialise the storage for the factors
      this->mFactors.resize(nFactors);

      // Compute the factors
      SplittingTools::factorizeNCpu(this->mNCpuFactors, nFactors, this->nCpu());

      // Filter the factors
      SplittingTools::filterFactors(this->mNCpuFactors, nFactors, this->nCpu(), true);
   }

   bool SplittingAlgorithm::useNextFactors()
   {
      // Get iterator through known factors
      auto it = this->mNCpuFactors.begin();

      // Check if list is empty
      if(it == this->mNCpuFactors.end())
      {
         return false;

      // Setup system with next factors
      } else
      {
         // Extract the next factors to try and remove them from list
         for(int i = 0; i < this->mFactors.size(); i++)
         {
            this->mFactors(i) = *it;
            it = this->mNCpuFactors.erase(it);
         }

         return true;
      }
   }

   double SplittingAlgorithm::communicationScore(SharedResolution spRes, ArrayI& details)
   {
      // The worst possible value is obtained for an all-to-all communication
      // at each (possible) communication step
      int worst = (this->dims()-1)*this->nCpu();

      // Initialise current structure score
      details.resize(spRes->cpu(0)->nDim()-1);
      details.setConstant(worst);

      if(spRes->cpu()->nDim()-1 == this->mFactors.size())
      {
         details = this->mFactors;
      }

      // Return ratio of both structures (higher is better)
      return static_cast<double>(worst)/static_cast<double>(this->mFactors.sum());
   }

   double SplittingAlgorithm::balancingScore(SharedResolution spRes, Array& balance)
   {
      // Storage for the per CPU loads for each dimension
      std::vector<std::map<int, double> >   loads;

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw std::logic_error("Requested computation of load balancing score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Loop over dimensions
         for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
         {
            // Create storage
            loads.push_back(std::map<int, double>());

            // Loop over CPUs
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // Initialise CPU load to zero
               loads.at(d)[cpu] = 0.0;

               // Loop over second dimension
               for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT2D>(); j++)
               {
                  // Increment load by 1
                  loads.at(d).find(cpu)->second += 1.0;
               }
            }
         }

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Loop over dimensions
         for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
         {
            // Create storage
            loads.push_back(std::map<int, double>());

            // Loop over CPUs
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // Initialise CPU fload to zero
               loads.at(d)[cpu] = 0.0;

               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT3D>(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT2D>(i); j++)
                  {
                     // Increment load by 1
                     loads.at(d).find(cpu)->second += 1.0;
                  }
               }
            }
         }
      }

      // Get total load
      double optimal = 0.0;
      Array perCpu(spRes->nCpu());

      // Loop over dimensions
      std::map<int, double>::const_iterator  it;
      for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
      {
         // Reset loads
         optimal = 0.0;
         perCpu.setConstant(0.0);

         for(it = loads.at(d).begin(); it != loads.at(d).end(); it++)
         {
            perCpu(it->first) += it->second;
            optimal += it->second;
         }

         // Convert total load to optimal load per CPU
         optimal = optimal/spRes->nCpu();

         // Load balance
         if(perCpu.minCoeff() > optimal)
         {
            balance(d) *= optimal/perCpu.maxCoeff();
         } else if(perCpu.maxCoeff() < optimal)
         {
            balance(d) *= perCpu.minCoeff()/optimal;
         } else
         {
            balance(d) *= std::min(perCpu.minCoeff()/optimal, optimal/perCpu.maxCoeff());
         }
      }

      // Compute score
      double score = 1.0;

      for(int i = 0; i < balance.size(); i++)
      {
         score *= balance(i);
      }

      return score;
   }

   void SplittingAlgorithm::buildCommunicationStructure(const int localId, SharedResolution spRes, std::map<Dimensions::Transform::Id,std::multimap<int,int> >& commStructure)
   {
      Profiler::RegionFixture<4> fix("Framework::LoadSplitter::SplittingAlgorithm::buildCommunicationStructure");

      // Clear the communication structure
      std::map<Dimensions::Transform::Id,std::multimap<int,int> >().swap(commStructure);

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw std::logic_error("Requested computation of communication structure score for 1D resolution!");
      }
      // Handle 2D resolution
      else if(spRes->cpu(0)->nDim() == 2)
      {
         details::buildCommunicationStructure2D(localId, spRes, commStructure);
      }
      // Handle 3D resolution
      else if(spRes->cpu(0)->nDim() == 3)
      {
         details::buildCommunicationStructure3D(localId, spRes, commStructure);
      }

      // Synchronize
      QuICCEnv().synchronize();
   }

}
}
