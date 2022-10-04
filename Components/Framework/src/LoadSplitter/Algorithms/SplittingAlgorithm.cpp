/** 
 * @file SplittingAlgorithm.cpp
 * @brief Source of the base of the implementation of the load splitting algorithms
 */

// Configuration includes
//

// System includes
//
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include <stdexcept>

#ifdef QUICC_MPI
#include <mpi.h>
#endif

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
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
      SplittingTools::filterFactors(this->mNCpuFactors, nFactors, this->nCpu());
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
      // Clear the communication structure
      std::map<Dimensions::Transform::Id,std::multimap<int,int> >().swap(commStructure);

      Dimensions::Transform::Id dimId;
      int i_;
      int j_;
      int k_;

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw std::logic_error("Requested computation of communication structure score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Simplify syntax
         typedef std::pair<int,int>   Coordinate;

         // Extract communication structure from resolution object
         std::set<Coordinate> bwdMap;
         std::set<Coordinate> fwdMap;

         // Storage for a coordinate
         Coordinate point;

         // Position iterator for insert calls
         std::set<Coordinate>::iterator   mapPos;

         // Loop over possible data exchanges 
         std::vector<Dimensions::Transform::Id> exchanges = {Dimensions::Transform::TRA1D};
         for(auto exId: exchanges)
         {
            // Create storage for structure
            commStructure.emplace(exId,std::multimap<int,int>());

            // initialise the position hint for inserts
            mapPos = bwdMap.begin();

            dimId = Dimensions::jump(exId,1);
            const auto& bwdTRes = *spRes->cpu()->dim(dimId);
            // Loop over second dimension
            for(int j = 0; j < bwdTRes.dim<Dimensions::Data::DAT2D>(); j++)
            {
               j_ = bwdTRes.idx<Dimensions::Data::DAT2D>(j);

               // Loop over backward dimension
               for(int i = 0; i < bwdTRes.dim<Dimensions::Data::DATB1D>(); i++)
               {
                  i_ = bwdTRes.idx<Dimensions::Data::DATB1D>(i);

                  // Generate point information
                  point = spRes->counter().makeKey(dimId, i_, j_);

                  // Get insertion position to use as next starting point to speed up insertion
                  mapPos = bwdMap.insert(mapPos,point);
               }
            }

            // Loop over CPUs
            MatrixI  matRemote;
            int matched = 0;
            int toMatch = -1;
            std::set<std::pair<int,int> > filter;
            dimId = exId;
            const auto& fwdTRes = *spRes->cpu()->dim(dimId);
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               matched = 0;

               // Local CPU
               if(cpu == localId)
               {
                  // Loop over second dimension
                  for(int j = 0; j < fwdTRes.dim<Dimensions::Data::DAT2D>(); j++)
                  {
                     j_ = fwdTRes.idx<Dimensions::Data::DAT2D>(j);

                     // Loop over forward dimension
                     for(int i = 0; i < fwdTRes.dim<Dimensions::Data::DATF1D>(); i++)
                     {
                        i_ = fwdTRes.idx<Dimensions::Data::DATF1D>(i);

                        // Generate point information
                        point = spRes->counter().makeKey(dimId, i_, j_);

                        // Look for same key in backward list
                        mapPos = bwdMap.find(point);

                        // Key was present, drop entry and extend filter
                        if(mapPos != bwdMap.end())
                        {
                           // Add corresponding communication edge to filter
                           filter.insert(std::make_pair(cpu, localId));

                           // Delete found coordinate
                           bwdMap.erase(mapPos);
                        } else
                        {
                           fwdMap.insert(point);
                        }
                     }
                  }

                  // Store size of forward coordinates
                  toMatch = fwdMap.size();

               #ifdef QUICC_MPI
                  // Convert coordinates set to matrix to send through MPI
                  matRemote.resize(2, fwdMap.size());
                  int i =0;
                  for(auto it = fwdMap.begin(); it != fwdMap.end(); ++it)
                  {
                     matRemote(0,i) = it->first;
                     matRemote(1,i) = it->second;
                     i++;
                  }

                  // Broadcast size
                  QuICCEnv().synchronize();
                  int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 711);

                  // Broadcast data
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                  QuICCEnv().check(ierr, 712);

               // Remote CPU   
               } else
               {
                  // Get size
                  QuICCEnv().synchronize();
                  int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 713);

                  // Get remote keys as matrix
                  matRemote.resize(2, toMatch);
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                  QuICCEnv().check(ierr, 714);

                  // Compare received data to stored indexes
                  for(int i = 0; i < toMatch; i++)
                  {
                     point = std::make_pair(matRemote(0,i), matRemote(1,i));

                     mapPos = bwdMap.find(point);

                     // Check if point is in backward map
                     if(mapPos != bwdMap.end())
                     {
                        // Add corresponding communication edge to filter
                        filter.insert(std::make_pair(cpu, localId));

                        // Delete found entry
                        bwdMap.erase(mapPos);

                        // Increase matched counter
                        matched++;
                     }
                  }
               }

               #else
               }
               #endif // QUICC_MPI
            }

            #ifdef QUICC_MPI
               // Gather total number of match entries
               QuICCEnv().synchronize();
               int ierr = MPI_Allreduce(MPI_IN_PLACE, &matched, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               QuICCEnv().check(ierr, 715);
            #endif // QUICC_MPI

            // Check that everything matched
            if(toMatch != matched)
            {
               throw std::logic_error("The computed index sets don't match!");
            }

            #ifdef QUICC_MPI
               // Store current filter
               MatrixI locFilter(2, filter.size());
               int i = 0;
               for(auto it = filter.begin(); it != filter.end(); ++it)
               {
                  locFilter(0, i) = it->first;
                  locFilter(1, i) = it->second;
                  i++;
               }

               // Gather full communication structure
               for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
               {
                  int filterSize = 0;
                  if(cpu == localId)
                  {
                     filterSize = locFilter.cols();

                     // Get size
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 716);

                     // Get remote keys as matrix
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(locFilter.data(), locFilter.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                     QuICCEnv().check(ierr, 717);

                  } else
                  {
                     // Get size
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 718);

                     // Get remote keys as matrix
                     matRemote.resize(2, filterSize);
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                     QuICCEnv().check(ierr, 719);

                     for(int i = 0; i < filterSize; ++i)
                     {
                        filter.insert(std::make_pair(matRemote(0,i), matRemote(1,i)));
                     }
                  }
               }
            #endif // QUICC_MPI

            // Store obtained minimized structure
            for(auto filIt = filter.begin(); filIt != filter.end(); filIt++)
            {
               commStructure.at(exId).insert(*filIt);
            }

            // Clear all the data for next loop
            bwdMap.clear();
            fwdMap.clear();
         }

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Simplify syntax
         typedef std::tuple<int,int,int>   Coordinate;

         // Extract communication structure from resolution object
         std::set<Coordinate> bwdMap;
         std::set<Coordinate> fwdMap;

         // Storage for a coordinate
         Coordinate point;

         // Position iterator for insert calls
         std::set<Coordinate>::iterator   mapPos;

         // Loop over possible data exchanges 
         std::vector<Dimensions::Transform::Id> exchanges = {Dimensions::Transform::TRA1D,Dimensions::Transform::TRA2D,Dimensions::Transform::SPECTRAL};
         for(auto exId: exchanges)
         {
            // Create storage for structure
            commStructure.emplace(exId, std::multimap<int,int>());

            // initialise the position hint for inserts
            mapPos = bwdMap.begin();

            dimId = Dimensions::jump(exId, 1);
            const auto& fwdTRes = *spRes->cpu()->dim(dimId);
            // Loop over third dimension
            for(int k = 0; k < fwdTRes.dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = fwdTRes.idx<Dimensions::Data::DAT3D>(k);

               // Loop over second dimension
               for(int j = 0; j < fwdTRes.dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = fwdTRes.idx<Dimensions::Data::DAT2D>(j,k);

                  // Loop over backward dimension
                  for(int i = 0; i < fwdTRes.dim<Dimensions::Data::DATB1D>(k); i++)
                  {
                     i_ = fwdTRes.idx<Dimensions::Data::DATB1D>(i,k);

                     // Generate point information
                     point = spRes->counter().makeKey(dimId, i_, j_, k_);

                     // Get insertion position to use as next starting point to speed up insertion
                     mapPos = bwdMap.insert(mapPos,point);
                  }
               }
            }

            // Loop over CPUs
            MatrixI  matRemote;
            int matched = 0;
            int toMatch = -1;
            std::set<std::pair<int,int> > filter;
            dimId = exId;
            const auto& bwdTRes = *spRes->cpu()->dim(dimId);
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               matched = 0;

               // Local CPU
               if(cpu == localId)
               {
                  // Loop over third dimension
                  for(int k = 0; k < bwdTRes.dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     k_ = bwdTRes.idx<Dimensions::Data::DAT3D>(k);

                     // Loop over second dimension
                     for(int j = 0; j < bwdTRes.dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        j_ = bwdTRes.idx<Dimensions::Data::DAT2D>(j,k);

                        // Loop over forward dimension
                        for(int i = 0; i < bwdTRes.dim<Dimensions::Data::DATF1D>(k); i++)
                        {
                           i_ = bwdTRes.idx<Dimensions::Data::DATF1D>(i,k);

                           // Generate point information
                           point = spRes->counter().makeKey(dimId, i_, j_, k_);

                           // Look for same key in backward list
                           mapPos = bwdMap.find(point);

                           // Key was present, drop entry and extend filter
                           if(mapPos != bwdMap.end())
                           {
                              // Add corresponding communication edge to filter
                              filter.insert(std::make_pair(cpu, localId));

                              // Delete found coordinate
                              bwdMap.erase(mapPos);
                           } else
                           {
                              fwdMap.insert(point);
                           }
                        }
                     }
                  }

                  // Store size of forward coordinates
                  toMatch = fwdMap.size();

               #ifdef QUICC_MPI
                  // Convert coordinates set to matrix to send through MPI
                  matRemote.resize(3, fwdMap.size());
                  int i =0;
                  for(auto it = fwdMap.begin(); it != fwdMap.end(); ++it)
                  {
                     matRemote(0,i) = std::get<0>(*it);
                     matRemote(1,i) = std::get<1>(*it);
                     matRemote(2,i) = std::get<2>(*it);
                     i++;
                  }

                  // Broadcast size
                  QuICCEnv().synchronize();
                  int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 720);

                  // Broadcast data
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                  QuICCEnv().check(ierr, 721);

               // Remote CPU   
               } else
               {
                  // Get size
                  QuICCEnv().synchronize();
                  int ierr = MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                  QuICCEnv().check(ierr, 722);

                  // Get remote keys as matrix
                  matRemote.resize(3, toMatch);
                  QuICCEnv().synchronize();
                  ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                  QuICCEnv().check(ierr, 723);

                  // Compare received data to stored indexes
                  for(int i = 0; i < toMatch; i++)
                  {
                     point = std::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));

                     mapPos = bwdMap.find(point);

                     // Check if point is in backward map
                     if(mapPos != bwdMap.end())
                     {
                        // Add corresponding communication edge to filter
                        filter.insert(std::make_pair(cpu, localId));

                        // Delete found entry
                        bwdMap.erase(mapPos);

                        // Increase matched counter
                        matched++;
                     }
                  }
               }

               #else
               }
               #endif // QUICC_MPI
            }

            #ifdef QUICC_MPI
               // Gather total number of match entries
               QuICCEnv().synchronize();
               int ierr = MPI_Allreduce(MPI_IN_PLACE, &matched, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               QuICCEnv().check(ierr, 724);
            #endif // QUICC_MPI

            // Check that everything matched
            if(toMatch != matched)
            {
               throw std::logic_error("The computed index sets don't match!");
            }

            #ifdef QUICC_MPI
               // Store current filter
               MatrixI locFilter(2, filter.size());
               int i = 0;
               for(auto it = filter.begin(); it != filter.end(); ++it)
               {
                  locFilter(0, i) = it->first;
                  locFilter(1, i) = it->second;
                  i++;
               }

               // Gather full communication structure
               for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
               {
                  int filterSize = 0;
                  if(cpu == localId)
                  {
                     filterSize = locFilter.cols();

                     // Get size
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 725);

                     // Get remote keys as matrix
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(locFilter.data(), locFilter.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                     QuICCEnv().check(ierr, 726);

                  } else
                  {
                     // Get size
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);
                     QuICCEnv().check(ierr, 727);

                     // Get remote keys as matrix
                     matRemote.resize(2, filterSize);
                     QuICCEnv().synchronize();
                     ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpu, MPI_COMM_WORLD); 
                     QuICCEnv().check(ierr, 728);

                     for(int i = 0; i < filterSize; ++i)
                     {
                        filter.insert(std::make_pair(matRemote(0,i), matRemote(1,i)));
                     }
                  }
               }
            #endif // QUICC_MPI

            // Store obtained minimized structure
            for(auto&& fId: filter)
            {
               commStructure.at(exId).insert(fId);
            }

            // Clear all the data for next loop
            bwdMap.clear();
            fwdMap.clear();
         }
      }

      // Synchronize
      QuICCEnv().synchronize();
   }

}
}
