/** 
 * @file RegularSHmIndexCounter.hpp
 * @brief Implementation of spherical harmonic index counter with m spectral and transform ordering with regular radial truncation
 */

#ifndef QUICC_REGULARSHMINDEXCOUNTER_HPP
#define QUICC_REGULARSHMINDEXCOUNTER_HPP

// Configuration includes
//

// System includes
//
#include <memory>
#include <tuple>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

namespace QuICC {

   /**
    * @brief Implementation of spherical harmonic index counter with m spectral ordering
    */ 
   class RegularSHmIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         RegularSHmIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~RegularSHmIndexCounter();

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version uses the internally stored simulation resolution
          *
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version reorders the input dimensions
          *
          * @param dims    Array of dimensions to reorder (1D, 2D, 3D, ...)
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Comput the offsets for the local modes
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Compute the offsets for the local modes by comparing to a reference simulation
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const;

         /**
          * @brief Generate index key as vector
          */
         virtual std::vector<int> makeVKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const;

         /**
          * @brief Generate index key
          */
         virtual std::tuple<int,int,int> makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const;
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a RegularSHmIndexCounter
   typedef std::shared_ptr<RegularSHmIndexCounter>   SharedSHmIndexCounter;

}

#endif // QUICC_REGULARSHMINDEXCOUNTER_HPP
