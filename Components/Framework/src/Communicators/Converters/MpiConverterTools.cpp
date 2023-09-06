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
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"
#include "QuICC/QuICCEnv.hpp"

namespace QuICC {

namespace Parallel {

   //
   // Three dimensional
   //

   void MpiConverterTools::buildLocalFwdMap(CoordinateMap&  rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim)
   {
      Parallel::PassthroughIndexConv idxConv;

      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildLocalMap3D<Dimensions::Data::DATF1D>(rLocalIdxMap, spRes, fwdDim, idxConv);
            break;
         case 2:
            MpiConverterTools::buildLocalMap2D<Dimensions::Data::DATF1D>(rLocalIdxMap, spRes, fwdDim, idxConv);
            break;
         case 1:
            MpiConverterTools::buildLocalMap1D<Dimensions::Data::DATF1D>(rLocalIdxMap, spRes, fwdDim, idxConv);
            break;
         default:
            throw std::logic_error("Tried to build local forward coordinate map for unimplemented dimension!");
      }
   }

   void MpiConverterTools::buildLocalBwdMap(CoordinateMap& rLocalIdxMap, SharedResolution spRes, const Dimensions::Transform::Id fwdDim, const IIndexConv& idxConv)
   {
      const Dimensions::Transform::Id tId = Dimensions::jump(fwdDim,1);

      switch(spRes->sim().ss().dimension())
      {
         case 3:
            MpiConverterTools::buildLocalMap3D<Dimensions::Data::DATB1D>(rLocalIdxMap, spRes, tId, idxConv);
            break;
         case 2:
            MpiConverterTools::buildLocalMap2D<Dimensions::Data::DATB1D>(rLocalIdxMap, spRes, tId, idxConv);
            break;
         case 1:
            MpiConverterTools::buildLocalMap1D<Dimensions::Data::DATB1D>(rLocalIdxMap, spRes, tId, idxConv);
            break;
         default:
            throw std::logic_error("Tried to build local forward coordinate map for unimplemented dimension!");
      }
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

}
}
