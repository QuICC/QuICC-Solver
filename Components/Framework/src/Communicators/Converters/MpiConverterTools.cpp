/**
 * @file MpiConverterTools.cpp
 * @brief Source of the tools for the MPI converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/MpiConverterTools.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Framework/MpiFramework.hpp"

namespace QuICC {

namespace Parallel {

   //
   // Three dimensional
   //

   void MpiConverterTools::buildLocalFwdMap(CoordinateMap&  rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildLocalFwdMap3D(rLocalIdxMap, spRes, fwdDim);
            break;
         case 2:
            MpiConverterTools::buildLocalFwdMap2D(rLocalIdxMap, spRes, fwdDim);
            break;
         case 1:
            MpiConverterTools::buildLocalFwdMap1D(rLocalIdxMap, spRes, fwdDim);
            break;
         default:
            throw std::logic_error("Tried to build local forward coordinate map for unimplemented dimension!");
      }
   }

   void MpiConverterTools::buildLocalBwdMap(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const IIndexConv& idxConv)
   {
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildLocalBwdMap3D(rLocalIdxMap, spRes, fwdDim, idxConv);
            break;
         case 2:
            MpiConverterTools::buildLocalBwdMap2D(rLocalIdxMap, spRes, fwdDim);
            break;
         case 1:
            MpiConverterTools::buildLocalBwdMap1D(rLocalIdxMap, spRes, fwdDim);
            break;
         default:
            throw std::logic_error("Tried to build local forward coordinate map for unimplemented dimension!");
      }
   }

   void MpiConverterTools::buildLocalFwdMap3D(CoordinateMap&  rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int k = 0; k < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over forward data dimension
            for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

               // Combine array indexes into coordinate tuple
               Coordinate coord = {i, j, k};

               // Create key as (1D, 2D, 3D)
               Coordinate key = spRes->counter().makeVKey(fwdDim, i_, j_, k_);

               // add key->coordinate to map
               rLocalIdxMap.insert(std::make_pair(key, coord));
            }
         }
      }
   }

   void MpiConverterTools::buildFwdCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes in next transform
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == MpiFramework::transformId(fwdDim))
      {
         // Loop over slow data dimension
         for(int k = 0; k < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // Extract "physical" index of slow data dimension
            k_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

            // Loop over middle data dimension
            for(int j = 0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               // Extract "physical" index of middle data dimension
               j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

               // Loop over backward data dimension
               for(int i = 0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
               {
                  // Extract "physical" index of backward data dimension
                  i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

                  // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
                  Coordinate key = spRes->counter().makeVKey(Dimensions::jump(fwdDim,1), i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send througg MPI
         matRemote.resize(3, remoteKeys.size());
         int col = 0;
         for(auto it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->at(0);
            matRemote(1,col) = it->at(1);
            matRemote(2,col) = it->at(2);
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 951);

         // Broadcast data
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 952);

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 953);

         // Resize storage
         matRemote.resize(3, nCoords);

         // Get remote keys as matrix
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 954);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            Coordinate key = {matRemote(0,i), matRemote(1,i), matRemote(2,i)};
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools::buildBwdCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == MpiFramework::transformId(fwdDim))
      {
         // Loop over slow data dimension
         for(int k = 0; k < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // Extract "physical" index of slow data dimension
            k_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

            // Loop over middle data dimension
            for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               // Extract "physical" index of middle data dimension
               j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

               // Loop over forward data dimension
               for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
               {
                  // Extract "physical" index of forward data dimension
                  i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

                  // Create key as (1D, 2D, 3D)
                  Coordinate key = spRes->counter().makeVKey(fwdDim, i_, j_, k_);

                  // Add key to remote set
                  remoteKeys.insert(key);
               }
            }
         }

         // Convert remote keys to matrix to send through MPI
         matRemote.resize(3, remoteKeys.size());
         int col = 0;
         for(auto it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->at(0);
            matRemote(1,col) = it->at(1);
            matRemote(2,col) = it->at(2);
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 961);

         // Broadcast data
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 962);

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 963);

         matRemote.resize(3, nCoords);

         // Get remot ekeys as matrix
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 964);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            Coordinate key = {matRemote(0,i), matRemote(1,i), matRemote(2,i)};
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools::extractShared(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, const std::set<Coordinate>& remoteKeys)
   {
      // List of local index keys
      std::set<Coordinate>  localKeys;

      // Extract the set of local keys
      for(auto mapIt = localIdxMap.begin(); mapIt != localIdxMap.end(); ++mapIt)
      {
         localKeys.insert(mapIt->first);
      }

      // Storage for the shared keys
      std::set<Coordinate> sharedKeys;

      // Create the list of common indexes
      std::set_intersection(localKeys.begin(), localKeys.end(), remoteKeys.begin(), remoteKeys.end(), std::inserter(sharedKeys, sharedKeys.begin()));

      // Clear the shared map
      sharedMap.clear();

      // Fill shared map
      for(auto sit = sharedKeys.begin(); sit != sharedKeys.end(); sit++)
      {
         sharedMap.insert(std::make_pair(*sit, localIdxMap.find(*sit)->second));
      }
   }

   void MpiConverterTools::buildLocalBwdMap3D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const IIndexConv& idxConv)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      //
      // Create the list of local indexes
      //

      // Loop over slow data dimension
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Extract "physical" index of slow data dimension
         k_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over backward data dimension
            for(int i = 0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

               // Combine array indexes into coordinate tuple
               Coordinate coord = {i + idxConv.centralPadding(i_, k), j, k};

               // Create key as (2D, 3D, 1D)
               Coordinate key = spRes->counter().makeVKey(Dimensions::jump(fwdDim,1), i_, j_, k_);

               // add key->coordinate to map
               rLocalIdxMap.insert(std::make_pair(key, coord));
            }
         }
      }
   }

   //
   // Two dimensional
   //

   void MpiConverterTools::buildLocalFwdMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_;

      //
      // Create the list of local indexes
      //

      // Loop over middle data dimension
      for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of middle data dimension
         j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

         // Loop over forward data dimension
         for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
         {
            // Extract "physical" index of forward data dimension
            i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i);

            // Combine array indexes into coordinate tuple
            Coordinate coord = {i, j};

            // Create key as (1D, 2D)
            Coordinate key = spRes->counter().makeVKey(fwdDim, i_, j_);

            // add key->coordinate to map
            rLocalIdxMap.insert(std::make_pair(key, coord));
         }
      }
   }

   void MpiConverterTools::buildFwdCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes in next transform
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == MpiFramework::transformId(fwdDim))
      {
         // Loop over middle data dimension
         for(int j=0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j);

            // Loop over backward data dimension
            for(int i=0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

               // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
               Coordinate key = spRes->counter().makeVKey(Dimensions::jump(fwdDim,1), i_, j_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

         // Convert remote keys to matrix to send througg MPI
         matRemote.resize(2, remoteKeys.size());
         int col = 0;
         for(auto it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->at(0);
            matRemote(1,col) = it->at(1);
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 951);

         // Broadcast data
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 952);

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 953);

         matRemote.resize(2, nCoords);

         // Get remote keys as matrix
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 954);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            Coordinate key = {matRemote(0,i), matRemote(1,i)};
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools::buildLocalBwdMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_;

      //
      // Create the list of local indexes
      //

      // Loop over middle data dimension
      for(int j = 0; j < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of middle data dimension
         j_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j);

         // Loop over backward data dimension
         for(int i = 0; i < spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
         {
            // Extract "physical" index of backward data dimension
            i_ = spRes->cpu()->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

            // Combine array indexes into coordinate tuple
            Coordinate coord = {i, j};

            // Create key as (2D, 3D, 1D)
            Coordinate key = spRes->counter().makeVKey(Dimensions::jump(fwdDim,1), i_, j_);

            // add key->coordinate to map
            rLocalIdxMap.insert(std::make_pair(key, coord));
         }
      }
   }

   void MpiConverterTools::buildBwdCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_, j_;

      // MPI error code
      int ierr;

      // Number of coordinates
      int nCoords = -1;

      //
      // Create the list of remote indexes
      //

      // Remote is also local CPU
      MatrixI  matRemote;
      if(cpuId == MpiFramework::transformId(fwdDim))
      {
         // Loop over middle data dimension
         for(int j = 0; j < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j);

            // Loop over forward data dimension
            for(int i = 0; i < spRes->cpu()->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
            {
               // Extract "physical" index of forward data dimension
               i_ = spRes->cpu()->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i);

               // Create key as (1D, 2D, 3D)
               Coordinate key = spRes->counter().makeVKey(fwdDim, i_, j_);

               // Add key to remote set
               remoteKeys.insert(key);
            }
         }

         // Convert remote keys to matrix to send through MPI
         matRemote.resize(2, remoteKeys.size());
         int col = 0;
         for(auto it = remoteKeys.begin(); it != remoteKeys.end(); ++it)
         {
            matRemote(0,col) = it->at(0);
            matRemote(1,col) = it->at(1);
            col++;
         }

         // Broadcast size
         nCoords = remoteKeys.size();
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 961);

         // Broadcast data
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 962);

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 963);

         matRemote.resize(2, nCoords);

         // Get remot ekeys as matrix
         MpiFramework::syncTransform(fwdDim);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, MpiFramework::transformComm(fwdDim));
         QuICCEnv().check(ierr, 964);

         // Convert matrix to remoteKeys set
         for(int i = 0; i < matRemote.cols(); i++)
         {
            Coordinate key = {matRemote(0,i), matRemote(1,i)};
            remoteKeys.insert(key);
         }
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   //
   // One dimensional
   //

   void MpiConverterTools::buildLocalFwdMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of local indexes
      //

      // Loop over forward data dimension
      for(int i=0; i < spRes->cpu(QuICCEnv().id())->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(QuICCEnv().id())->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i);

         // Set coordindate
         Coordinate coord = {i};

         // Set key
         Coordinate key = {i_};

         // add key->coordinate to map
         rLocalIdxMap.insert(std::make_pair(key, coord));
      }
   }

   void MpiConverterTools::buildFwdCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of remote indexes in next transform
      //

      // Loop over backward data dimension
      for(int i=0; i < spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

         // Set key
         Coordinate key = {i_};

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

   void MpiConverterTools::buildLocalBwdMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of local indexes
      //

      // Loop ver backward data dimension
      for(int i=0; i < spRes->cpu(QuICCEnv().id())->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(QuICCEnv().id())->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i);

         // Set coordinate
         Coordinate coord = {i};

         // Set key
         Coordinate key = {i_};

         // add key->coordinate to map
         rLocalIdxMap.insert(std::make_pair(key, coord));
      }

   }

   void MpiConverterTools::buildBwdCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of remote indexes
      //

      // Loop over forward data dimension
      for(int i=0; i < spRes->cpu(cpuId)->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(); ++i)
      {
         // Extract "physical" index of forward data dimension
         i_ = spRes->cpu(cpuId)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATF1D>(i);

         // Set key
         Coordinate key = {i_};

         // Add key to remote set
         remoteKeys.insert(key);
      }

      // Extract map of shared indexes (stored as keys)
      sharedMap.clear();
      MpiConverterTools::extractShared(sharedMap, localIdxMap, remoteKeys);
   }

}
}
