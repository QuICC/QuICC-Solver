/** 
 * @file DoublePeriodicIndexCounter.hpp
 * @brief Implementation of index counter for doubly periodic schemes
 */

#ifndef QUICC_DOUBLEPERIODICINDEXCOUNTER_HPP
#define QUICC_DOUBLEPERIODICINDEXCOUNTER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

namespace QuICC {

   /**
    * @brief Implementation of regular index counter
    */ 
   class DoublePeriodicIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         DoublePeriodicIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~DoublePeriodicIndexCounter() = default;

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
          * @brief Compute the offsets for the local modes
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Compute the offsets for the local modes by comparing to a reference simulation
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const;
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a DoublePeriodicIndexCounter
   typedef std::shared_ptr<DoublePeriodicIndexCounter>   SharedDoublePeriodicIndexCounter;

} // QuICC

#endif // QUICC_DOUBLEPERIODICINDEXCOUNTER_HPP
