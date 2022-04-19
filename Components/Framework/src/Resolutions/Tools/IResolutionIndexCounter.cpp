/** 
 * @file IResolutionIndexCounter.cpp
 * @brief Source of regular index counter
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Resolutions/Tools/IResolutionIndexCounter.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

   IResolutionIndexCounter::IResolutionIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IndexCounter(), mspSim(spSim), mspCpu(spCpu)
   {
   }

   IResolutionIndexCounter::~IResolutionIndexCounter()
   {
   }

   int IResolutionIndexCounter::dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat) const
   {
      return this->mspSim->dim(simId, spaceId);
   }

   ArrayI IResolutionIndexCounter::dimensions(const Dimensions::Space::Id spaceId, const MHDFloat idx) const
   {
      auto dim = this->mspSim->dimensions(spaceId);
      for(int i = 0; i < dim.size(); i++)
      {
         dim(i) = this->dim(static_cast<Dimensions::Simulation::Id>(i), spaceId, idx);
      }
      return dim;
   }
}
