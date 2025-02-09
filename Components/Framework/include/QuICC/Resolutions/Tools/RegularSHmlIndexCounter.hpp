/** 
 * @file RegularSHmlIndexCounter.hpp
 * @brief Implementation of spherical harmonic index counter with m spectral ordering, regular radial truncation, and l transform ordering
 */

#ifndef QUICC_REGULARSHMLINDEXCOUNTER_HPP
#define QUICC_REGULARSHMLINDEXCOUNTER_HPP

// System includes
//
#include <memory>
#include <tuple>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

namespace QuICC {

   /**
    * @brief Implementation of spherical harmonic index counter with m spectral ordering, regular radial truncation and l transform ordering
    */ 
   class RegularSHmlIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         RegularSHmlIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~RegularSHmlIndexCounter() = default;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version uses the internally stored simulation resolution
          *
          * @param spaceId Spacial the resolution represent
          */
         ArrayI orderedDimensions(const Dimensions::Space::Id spaceId) const final;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version reorders the input dimensions
          *
          * @param dims    Array of dimensions to reorder (1D, 2D, 3D, ...)
          * @param spaceId Spacial the resolution represent
          */
         ArrayI orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const final;

         /**
          * @brief Comput the offsets for the local modes
          */
         void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const final;

         /**
          * @brief Compute the offsets for the local modes by comparing to a reference simulation
          */
         void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const final;

         /**
          * @brief Generate index key as vector
          */
         std::vector<int> makeVKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const final;

         /**
          * @brief Generate index key
          */
         std::tuple<int,int,int> makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const final;
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a RegularSHmlIndexCounter
   typedef std::shared_ptr<RegularSHmlIndexCounter>   SharedSHmIndexCounter;

} // QuICC

#endif // QUICC_REGULARSHMLINDEXCOUNTER_HPP
