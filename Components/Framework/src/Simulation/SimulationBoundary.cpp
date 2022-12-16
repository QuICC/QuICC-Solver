/**
 * @file SimulationBoundary.cpp
 * @brief Implementation of a general simulation control structure
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Simulation/SimulationBoundary.hpp"

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"

namespace QuICC {

   SimulationBoundary::SimulationBoundary(const std::map<std::string,std::size_t>& bcIds)
   {
      this->convert(bcIds);
   }

   SimulationBoundary::~SimulationBoundary()
   {
   }

   void SimulationBoundary::convert(const std::map<std::string,std::size_t>& bcIds)
   {
      for(auto mapIt = bcIds.cbegin(); mapIt != bcIds.cend(); ++mapIt)
      {
         std::size_t pn = Hasher::makeId(mapIt->first);

         this->mBcs.insert(std::make_pair(pn, mapIt->second));
      }

   }

   const std::map<std::size_t,std::size_t>& SimulationBoundary::map() const
   {
      return this->mBcs;
   }

   std::size_t SimulationBoundary::bcId(const std::size_t id) const
   {
      assert(this->mBcs.count(id) > 0);
      return this->mBcs.at(id);
   }

}
