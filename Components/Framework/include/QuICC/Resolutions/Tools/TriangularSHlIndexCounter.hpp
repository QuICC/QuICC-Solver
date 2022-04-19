/** 
 * @file TriangularSHlIndexCounter.hpp
 * @brief Implementation of spherical harmonic index counter with l spectral ordering and triangular radial truncation
 */

#ifndef QUICC_TRIANGULARSHLINDEXCOUNTER_HPP
#define QUICC_TRIANGULARSHLINDEXCOUNTER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

namespace QuICC {

   /**
    * @brief Implementation of spherical harmonic index counter with l spectral ordering
    */ 
   class TriangularSHlIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         TriangularSHlIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~TriangularSHlIndexCounter();

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         virtual int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat idx) const;

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
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a TriangularSHlIndexCounter
   typedef std::shared_ptr<TriangularSHlIndexCounter>   SharedSHlIndexCounter;

}

#endif // QUICC_TRIANGULARSHLINDEXCOUNTER_HPP
