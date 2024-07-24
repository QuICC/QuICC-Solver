/**
 * @file MpiConverterTools.hpp
 * @brief Implementation of some tools used by the MPI converter
 */

#ifndef QUICC_PARALLEL_MPICONVERTERTOOLS_HPP
#define QUICC_PARALLEL_MPICONVERTERTOOLS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <cassert>
#include <set>
#include <map>
#include <tuple>
#include <memory>

// External includes
//

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Types/Typedefs.hpp"
#include "Environment/MpiTypes.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Communicators/Converters/IIndexConv.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Specialisation of the tools used by the MPI converter for 3D simulation.
    */
   class MpiConverterTools
   {
      public:
         /// Typedef for a n indexes specifing a point
         typedef std::vector<int> Coordinate;
         /// Typedef for vector of coordinates
         typedef std::vector<Coordinate> CoordinateVector;
         /// Typedef for map from between coordinates
         typedef std::map<Coordinate,Coordinate> CoordinateMap;

         /**
          * @brief Build local cpu map of forward data
          */
         static void buildLocalFwdMap(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build local cpu map of backward data
          */
         static void buildLocalBwdMap(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const IIndexConv& idxConv);

         /**
          * @brief Build a MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData, Dimensions::Data::Id TDataId = Dimensions::Data::DATF1D> static MPI_Datatype buildDatatype(const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const Dimensions::Transform::Id keyId, TData& data, const int cpuId);

         /**
          * @brief Build manual packing description for Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData, Dimensions::Data::Id TDataId = Dimensions::Data::DATF1D> static void buildDatatype(CoordinateVector& coords, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const Dimensions::Transform::Id keyId, TData& data, const int cpuId);

         /**
          * @brief Pack data into buffer
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static void pack(TData* rBuffer, int &position, const typename Framework::Selector::ScalarField<TData> &data, const CoordinateVector& coords);

         /**
          * @brief Pack data into buffer
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static void pack(char* rBuffer, int &position, const typename Framework::Selector::ScalarField<TData> &data, const CoordinateVector& coords);

         /**
          * @brief Unpack data from buffer
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static void unpack(typename Framework::Selector::ScalarField<TData> &rData, const CoordinateVector& coords, const char* buffer, int &position);

         /**
          * @brief Unpack data from buffer
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static void unpack(typename Framework::Selector::ScalarField<TData> &rData, const CoordinateVector& coords, const TData* buffer, int &position);

      protected:

      private:
         /**
          * @brief Build local cpu map of forward 3D data
          */
         template <Dimensions::Data::Id TDataId> static void buildLocalMap3D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv);

         /**
          * @brief Build local cpu map of 2D data
          */
         template <Dimensions::Data::Id TDataId> static void buildLocalMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv);
         /**
          * @brief Build local cpu map of 1D data
          */
         template <Dimensions::Data::Id TDataId> static void buildLocalMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv);

         /**
          * @brief Build a forward 3D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param tId     Dimension index for transform
          * @param keyId   Dimension index for data keys
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <Dimensions::Data::Id TDataId> static void buildCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId);

         /**
          * @brief Build a forward 2D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param tId     Dimension index for transform
          * @param keyId   Dimension index for data keys
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <Dimensions::Data::Id TDataId> static void buildCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId);

         /**
          * @brief Build a forward 1D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param tId     Dimension index for transform
          * @param keyId   Dimension index for data keys
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <Dimensions::Data::Id TDataId> static void buildCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId);

         /**
          * @brief Extract shared indexes
          *
          * @param sharedMap     Storage for the shared index map
          * @param localIdxMap   Local node key to indexes map
          * @param remoteKeys    Remote node index keys
          */
         static void extractShared(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, const std::set<Coordinate>& remoteKeys);

         /**
          * @brief Create type
          *
          * @param data       The concerned data
          * @param sharedMap  Shared index map
          * @param type       MPI datatype storage
          */
         template <typename TData> static MPI_Datatype buildType(typename Framework::Selector::ScalarField<TData> &data, const CoordinateMap& sharedMap);

         /**
          * @brief Constructor
          */
         MpiConverterTools();

         /**
          * @brief Destructor
          */
         ~MpiConverterTools();

   };

   template <typename TData>
      void MpiConverterTools::pack(char* rBuffer, int &position, const typename Framework::Selector::ScalarField<TData> &data, const CoordinateVector &coords)
   {
      // Create MPI displacement list
      for(auto it = coords.cbegin(); it != coords.cend(); ++it)
      {
         memcpy(rBuffer+position, &data.point(*it), sizeof(TData));

         // Increment datatype size
         position += sizeof(TData);
      }
   }

   template <typename TData>
      void MpiConverterTools::pack(TData* rBuffer, int &position, const typename Framework::Selector::ScalarField<TData> &data, const CoordinateVector &coords)
   {
      // Create MPI displacement list
      for(auto it = coords.cbegin(); it != coords.cend(); ++it)
      {
         *(rBuffer+position) = data.point(*it);

         // Increment datatype size
         position++;
      }
   }

   template <typename TData>
      void MpiConverterTools::unpack(typename Framework::Selector::ScalarField<TData> &rData, const CoordinateVector& coords, const char* buffer, int &position)
   {
      // Create MPI displacement list
      for(auto it = coords.cbegin(); it != coords.cend(); ++it)
      {
         memcpy(&rData.rPoint(*it), buffer+position, sizeof(TData));

         // Increment datatype size
         position += sizeof(TData);
      }
   }

   template <typename TData>
      void MpiConverterTools::unpack(typename Framework::Selector::ScalarField<TData> &rData, const CoordinateVector& coords, const TData* buffer, int &position)
   {
      // Create MPI displacement list
      for(auto it = coords.cbegin(); it != coords.cend(); ++it)
      {
         rData.rPoint(*it) = *(buffer + position);

         // Increment datatype size
         position++;
      }
   }

   template <typename TData>
      MPI_Datatype MpiConverterTools::buildType(typename Framework::Selector::ScalarField<TData> &data, const CoordinateMap& sharedMap)
   {
      // MPI error code
      int ierr;

      // Get the number of elements
      int nElements = sharedMap.size();

      // Prepare data required for MPI datatype
      MPI_Aint *displ = new MPI_Aint[nElements];
      int      *blocks = new int[nElements];

      // Prepare data required for MPI datatype
      MPI_Aint    base;
      MPI_Aint    element;

      // Set the base of the datatype to the start of the data matrix
      MPI_Get_address(data.rData().data(),&element);
      base = element;

      // Create MPI displacement list
      int tot = 0;
      for(auto it = sharedMap.cbegin(); it != sharedMap.cend(); ++it)
      {
         // Get address of stored coordinates
         MPI_Get_address(&data.rPoint(it->second), &element);

         // Fill datatype information
         displ[tot] = element - base;
         blocks[tot] = 1;

         // Increment datatype size
         tot++;
      }

      // Create MPI datatype
      MPI_Datatype type;
      //ierr = MPI_Type_create_hindexed(nElements, blocks, displ, MpiTypes::type<TData>(), &type);
      ierr = MPI_Type_create_hindexed_block(nElements, 1, displ, Environment::MpiTypes::type<TData>(), &type);
      QuICCEnv().check(ierr, 971);

      // Commit MPI datatype
      ierr = MPI_Type_commit(&type);
      QuICCEnv().check(ierr, 972);

      // Free allocated memory
      delete[] displ;
      delete[] blocks;

      return type;
   }

   template <typename TData, Dimensions::Data::Id TDataId>
      void MpiConverterTools::buildDatatype(CoordinateVector& coords, const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const Dimensions::Transform::Id keyId, TData& data, const int cpuId)
   {
      CoordinateMap sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildCpuMap3D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         case 2:
            MpiConverterTools::buildCpuMap2D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         case 1:
            MpiConverterTools::buildCpuMap1D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         default:
            throw std::logic_error("Tried to build CPU coordinate map for unimplemented dimension!");
      }

      coords.clear();
      coords.reserve(sharedMap.size());
      for(auto it = sharedMap.cbegin(); it != sharedMap.cend(); ++it)
      {
         coords.push_back(it->second);
      }
   }

   template <typename TData, Dimensions::Data::Id TDataId>
      MPI_Datatype MpiConverterTools::buildDatatype(const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const Dimensions::Transform::Id keyId, TData& data, const int cpuId)
   {
      CoordinateMap sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildCpuMap3D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         case 2:
            MpiConverterTools::buildCpuMap2D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         case 1:
            MpiConverterTools::buildCpuMap1D<TDataId>(sharedMap, localIdxMap, spRes, fwdDim, keyId, cpuId);
            break;
         default:
            throw std::logic_error("Tried to build CPU coordinate map for unimplemented dimension!");
      }

      // Create the datatype
      MPI_Datatype type = MpiConverterTools::buildType(data, sharedMap);

#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = 22*Debug::MemorySize<MPI_Aint>::BYTES*sharedMap.size();
      StorageProfilerMacro_update(StorageProfilerMacro::MPIBTYPES, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);
#endif //QUICC_STORAGEPROFILE

      return type;
   }

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildLocalMap3D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_, k_;

      //
      // Create the list of local indexes
      //

      const auto& tRes = *spRes->cpu()->dim(tId);
      // Loop over slow data dimension
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // Global index
         // Extract "physical" index of slow data dimension
         k_ = tRes.idx<Dimensions::Data::DAT3D>(k);

         // Loop over middle data dimension
         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);

            // Loop over fast data dimension
            for(int i = 0; i < tRes.dim<TDataId>(j,k); ++i)
            {
               // Extract "physical" index of fast data dimension
               i_ = tRes.idx<TDataId>(i,j,k);

               // Combine array indexes into coordinate tuple
               Coordinate coord = {i + idxConv.centralPadding(i_, k), j, k};

               // Create key
               Coordinate key = spRes->counter().makeVKey(tId, i_, j_, k_);

               // add key->coordinate to map
               rLocalIdxMap.insert(std::make_pair(key, coord));
            }
         }
      }
   }

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildLocalMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_, j_;

      //
      // Create the list of local indexes
      //

      const auto& tRes = *spRes->cpu()->dim(tId);
      // Loop over middle data dimension
      for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(); ++j)
      {
         // Extract "physical" index of middle data dimension
         j_ = tRes.idx<Dimensions::Data::DAT2D>(j);

         // Loop over fast data dimension
         for(int i = 0; i < tRes.dim<TDataId>(); ++i)
         {
            // Extract "physical" index of fast data dimension
            i_ = tRes.idx<TDataId>(i);

            // Combine array indexes into coordinate tuple
            Coordinate coord = {i, j};

            // Create key
            Coordinate key = spRes->counter().makeVKey(tId, i_, j_);

            // add key->coordinate to map
            rLocalIdxMap.insert(std::make_pair(key, coord));
         }
      }
   }

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildLocalMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const IIndexConv& idxConv)
   {
      // Make sure key map is empty
      rLocalIdxMap.clear();

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of local indexes
      //

      // Loop ver backward fast dimension
      for(int i=0; i < spRes->cpu(QuICCEnv().id())->dim(tId)->dim<TDataId>(); ++i)
      {
         // Extract "physical" index of fast data dimension
         i_ = spRes->cpu(QuICCEnv().id())->dim(tId)->idx<TDataId>(i);

         // Set coordinate
         Coordinate coord = {i};

         // Set key
         Coordinate key = {i_};

         // add key->coordinate to map
         rLocalIdxMap.insert(std::make_pair(key, coord));
      }

   }

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId)
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
      if(cpuId == QuICCEnv().id(tId))
      {
         const auto& tRes = *spRes->cpu()->dim(keyId);
         // Loop over slow data dimension
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // Extract "physical" index of slow data dimension
            k_ = tRes.idx<Dimensions::Data::DAT3D>(k);

            // Loop over middle data dimension
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               // Extract "physical" index of middle data dimension
               j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);

               // Loop over fast data dimension
               for(int i = 0; i < tRes.dim<TDataId>(j,k); ++i)
               {
                  // Extract "physical" index of fast data dimension
                  i_ = tRes.idx<TDataId>(i,j,k);

                  // Create key
                  Coordinate key = spRes->counter().makeVKey(keyId, i_, j_, k_);

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
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Broadcast of sizes in 3D remote index map failed");

         // Broadcast data
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Broadcast of 3D remote index map failed");

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Receiving broadcast of sizes of 3D remote index map failed");

         matRemote.resize(3, nCoords);

         // Get remote keys as matrix
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Receiving broadcast of 3D remote index map failed");

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

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId)
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
      if(cpuId == QuICCEnv().id(tId))
      {
         const auto& tRes = *spRes->cpu()->dim(keyId);
         // Loop over middle data dimension
         for(int j=0; j < tRes.dim<Dimensions::Data::DAT2D>(); ++j)
         {
            // Extract "physical" index of middle data dimension
            j_ = tRes.idx<Dimensions::Data::DAT2D>(j);

            // Loop over backward data dimension
            for(int i=0; i < tRes.dim<TDataId>(); ++i)
            {
               // Extract "physical" index of backward data dimension
               i_ = tRes.idx<TDataId>(i);

               // Create key as (2D, 3D, 1D) indexes (i.e. data gets transposed during communication)
               Coordinate key = spRes->counter().makeVKey(keyId, i_, j_);

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
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Broadcast of sizes in 2D remote index map failed");

         // Broadcast data
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Broadcast of 2D remote index map failed");

      // Remote CPU needs to generate list
      } else
      {
         // Get size
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(&nCoords, 1, MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Receiving broadcast of sizes of 2D remote index map failed");

         matRemote.resize(2, nCoords);

         // Get remote keys as matrix
         QuICCEnv().synchronize(tId);
         ierr = MPI_Bcast(matRemote.data(), matRemote.size(), MPI_INT, cpuId, QuICCEnv().comm(tId));
         QuICCEnv().check(ierr, "Receiving broadcast of 2D remote index map failed");

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

   template <Dimensions::Data::Id TDataId> void MpiConverterTools::buildCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id tId, const Dimensions::Transform::Id keyId, const int cpuId)
   {
      // List of remote keys
      std::set<Coordinate>  remoteKeys;

      // Storage for the simulation wide indexes
      int i_;

      //
      // Create the list of remote indexes in next transform
      //

      // Loop over backward data dimension
      for(int i=0; i < spRes->cpu(cpuId)->dim(keyId)->dim<TDataId>(); ++i)
      {
         // Extract "physical" index of backward data dimension
         i_ = spRes->cpu(cpuId)->dim(keyId)->idx<TDataId>(i);

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

#endif // QUICC_PARALLEL_MPICONVERTERTOOLS_HPP
