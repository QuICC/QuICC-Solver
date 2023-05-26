/** 
 * @file TrapezoidalSHlIndexCounter.hpp
 * @brief Implementation of spherical harmonic index counter with l spectral ordering and trapezoidal radial truncation
 */

#ifndef QUICC_TRAPEZOIDALSHLINDEXCOUNTER_HPP
#define QUICC_TRAPEZOIDALSHLINDEXCOUNTER_HPP

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
    * @brief Implementation of spherical harmonic index counter with l spectral ordering
    */ 
   class TrapezoidalSHlIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         TrapezoidalSHlIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~TrapezoidalSHlIndexCounter() = default;

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat idx) const final;

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
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a TrapezoidalSHlIndexCounter
   typedef std::shared_ptr<TrapezoidalSHlIndexCounter>   SharedTrapezoidalSHlIndexCounter;

} // QuICC

#endif // QUICC_TRAPEZOIDALSHLINDEXCOUNTER_HPP
