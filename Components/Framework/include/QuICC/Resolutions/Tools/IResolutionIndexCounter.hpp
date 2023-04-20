/** 
 * @file IResolutionIndexCounter.hpp
 * @brief Implementation of index counter with resolution objects
 */

#ifndef QUICC_IRESOLUTIONINDEXCOUNTER_HPP
#define QUICC_IRESOLUTIONINDEXCOUNTER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Resolutions/SimulationResolution.hpp"
#include "QuICC/Resolutions/CoreResolution.hpp"

namespace QuICC {

   /**
    * @brief Implementation ofr index counter with resolution objectbs
    */ 
   class IResolutionIndexCounter: public IndexCounter
   {
      public:
         /**
          * @brief Constructor
          */
         IResolutionIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu);

         /**
          * @brief Empty destructor
          */
         virtual ~IResolutionIndexCounter() = default;

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         virtual int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat idx) const;

         /**
          * @brief Get dimensions
          *
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         virtual ArrayI dimensions(const Dimensions::Space::Id spaceId, const MHDFloat idx) const;

      protected:
         /**
          * @brief Local copy of the simulation resolution
          */
         SharedCSimulationResolution mspSim;

         /**
          * @brief Local copy of the local resolution
          */
         SharedCCoreResolution mspCpu;

      private:

   };

} // QuICC

#endif // QUICC_IRESOLUTIONINDEXCOUNTER_HPP
