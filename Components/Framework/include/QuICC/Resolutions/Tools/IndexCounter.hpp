/** 
 * @file IndexCounter.hpp
 * @brief Implementation of base class for a generalized index counter
 */

#ifndef QUICC_INDEXCOUNTER_HPP
#define QUICC_INDEXCOUNTER_HPP

// System includes
//
#include <tuple>
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Io/Hdf5Typedefs.hpp"

namespace QuICC {

   // Forward declaration to avoid header loop
   class SimulationResolution;
   typedef std::shared_ptr<const SimulationResolution> SharedCSimulationResolution;

   /**
    * @brief Implementation of base class for a generalized index counter
    */ 
   class IndexCounter
   {
      public:
         /// Typedef for the offsets
         typedef Io::QuICC_hsize_t OffsetType;

         /**
          * @brief Constructor
          */
         IndexCounter();

         /**
          * @brief Empty destructor
          */
         virtual ~IndexCounter() = default;

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         virtual int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat idx) const = 0;

         /**
          * @brief Get dimensions
          *
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         virtual ArrayI dimensions(const Dimensions::Space::Id spaceId, const MHDFloat idx) const = 0;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version uses the internally stored simulation resolution
          *
          * @param spaceId Space the resolution represent
          */
         virtual ArrayI orderedDimensions(const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Reorder dimensions from fast to slow
          *
          * This version reorders the input dimensions
          *
          * @param dims    Array of dimensions to reorder (1D, 2D, 3D, ...)
          * @param spaceId Spacial the resolution represent
          */
         virtual ArrayI orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Compute the offset for local modes
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const = 0;

         /**
          * @brief Compute the offset for local modes by comparing to a reference simulation
          */
         virtual void computeOffsets(std::vector<OffsetType>& blocks, std::vector<std::vector<OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const = 0;

         /**
          * @brief Generate index key as vector
          */
         virtual std::vector<int> makeVKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const;

         /**
          * @brief Generate index key as vector
          */
         virtual std::vector<int> makeVKey(const Dimensions::Transform::Id id, const int i, const int j) const;

         /**
          * @brief Generate index key
          */
         virtual std::tuple<int,int,int> makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const;

         /**
          * @brief Generate index key
          */
         virtual std::pair<int,int> makeKey(const Dimensions::Transform::Id id, const int i, const int j) const;
         
      protected:

      private:

   };

   /// Typedef for an smart reference counting pointer for a IndexCounter
   typedef std::shared_ptr<IndexCounter>   SharedIndexCounter;

} // QuICC

#endif // QUICC_INDEXCOUNTER_HPP
