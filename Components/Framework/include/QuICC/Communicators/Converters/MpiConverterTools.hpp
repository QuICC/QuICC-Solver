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
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/MpiTypes.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Communicators/Converters/IIndexConv.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"
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
          * @brief Build a forward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildFwdDatatype(const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId);

         /**
          * @brief Build a backward MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static MPI_Datatype buildBwdDatatype(const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId);

         /**
          * @brief Build manual packing description for forward Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static void buildFwdDatatype(CoordinateVector& coords, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId);


         /**
          * @brief Build manual packing description for backward Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         template <typename TData> static void buildBwdDatatype(CoordinateVector& coords, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId);

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
         static void buildLocalFwdMap3D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build local cpu map of forward 2D data
          */
         static void buildLocalFwdMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build local cpu map of forward 1D data
          */
         static void buildLocalFwdMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build local cpu map of backward 3D data
          */
         static void buildLocalBwdMap3D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const IIndexConv& idxConv);

         /**
          * @brief Build local cpu map of backward 2D data
          */
         static void buildLocalBwdMap2D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build local cpu map of backward 1D data
          */
         static void buildLocalBwdMap1D(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim);

         /**
          * @brief Build a forward 3D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildFwdCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

         /**
          * @brief Build a forward 2D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildFwdCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

         /**
          * @brief Build a forward 1D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildFwdCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

         /**
          * @brief Build a backward 3D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildBwdCpuMap3D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

         /**
          * @brief Build a backward 2D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildBwdCpuMap2D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

         /**
          * @brief Build a backward 1D MPI Datatype
          *
          * @param spRes   Shared Resolution
          * @param fwdDim  Dimension index for forward transform
          * @param data    Input data
          * @param cpuId   ID of the CPU
          */
         static void buildBwdCpuMap1D(CoordinateMap& sharedMap, const CoordinateMap& localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const int cpuId);

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
      ierr = MPI_Type_create_hindexed_block(nElements, 1, displ, MpiTypes::type<TData>(), &type);
      QuICCEnv().check(ierr, 971);

      // Commit MPI datatype
      ierr = MPI_Type_commit(&type);
      QuICCEnv().check(ierr, 972);

      // Free allocated memory
      delete[] displ;
      delete[] blocks;

      return type;
   }

   template <typename TData>
      void MpiConverterTools::buildFwdDatatype(CoordinateVector& coords, const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId)
   {
      CoordinateMap  sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildFwdCpuMap3D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 2:
            MpiConverterTools::buildFwdCpuMap2D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 1:
            MpiConverterTools::buildFwdCpuMap1D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         default:
            throw std::logic_error("Tried to build CPU forward coordinate map for unimplemented dimension!");
      }

      coords.clear();
      coords.reserve(sharedMap.size());
      for(auto it = sharedMap.cbegin(); it != sharedMap.cend(); ++it)
      {
         coords.push_back(it->second);
      }
   }

   template <typename TData>
      void MpiConverterTools::buildBwdDatatype(CoordinateVector& coords, const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId)
   {
      CoordinateMap sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildBwdCpuMap3D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 2:
            MpiConverterTools::buildBwdCpuMap2D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 1:
            MpiConverterTools::buildBwdCpuMap1D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         default:
         throw std::logic_error("Tried to build CPU forward coordinate map for unimplemented dimension!");
      }

      coords.clear();
      coords.reserve(sharedMap.size());
      for(auto it = sharedMap.cbegin(); it != sharedMap.cend(); ++it)
      {
         coords.push_back(it->second);
      }
   }

   template <typename TData>
      MPI_Datatype MpiConverterTools::buildBwdDatatype(const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId)
   {
      CoordinateMap sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildBwdCpuMap3D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 2:
            MpiConverterTools::buildBwdCpuMap2D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 1:
            MpiConverterTools::buildBwdCpuMap1D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         default:
            throw std::logic_error("Tried to build CPU forward coordinate map for unimplemented dimension!");
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

   template <typename TData>
      MPI_Datatype MpiConverterTools::buildFwdDatatype(const CoordinateMap&  localIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, typename Framework::Selector::ScalarField<TData> &data, const int cpuId)
   {
      CoordinateMap  sharedMap;
      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildFwdCpuMap3D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 2:
            MpiConverterTools::buildFwdCpuMap2D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         case 1:
            MpiConverterTools::buildFwdCpuMap1D(sharedMap, localIdxMap, spRes, fwdDim, cpuId);
            break;
         default:
            throw std::logic_error("Tried to build CPU forward coordinate map for unimplemented dimension!");
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
}
}

#endif // QUICC_PARALLEL_MPICONVERTERTOOLS_HPP
