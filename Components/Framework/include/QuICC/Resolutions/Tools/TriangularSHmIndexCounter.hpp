/**
 * @file TriangularSHmIndexCounter.hpp
 * @brief Implementation of spherical harmonic index counter with m spectral ordering with triangular radial truncation
 */

#ifndef QUICC_TRIANGULARSHMINDEXCOUNTER_HPP
#define QUICC_TRIANGULARSHMINDEXCOUNTER_HPP

// System includes
//
#include <memory>
#include <tuple>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

namespace QuICC {

   /**
    * @brief Implementation of spherical harmonic index counter with m spectral ordering
    */
   class TriangularSHmIndexCounter: public IResolutionIndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         TriangularSHmIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         ~TriangularSHmIndexCounter() = default;

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

   /// Typedef for an smart reference counting pointer for a TriangularSHmIndexCounter
   typedef std::shared_ptr<TriangularSHmIndexCounter>   SharedSHmIndexCounter;

} // QuICC

#endif // QUICC_TRIANGULARSHMINDEXCOUNTER_HPP
