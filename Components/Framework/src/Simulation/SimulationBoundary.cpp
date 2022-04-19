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

   SimulationBoundary::SimulationBoundary(const std::map<std::string,int>& bcIds)
   {
      this->convert(bcIds);
   }

   SimulationBoundary::~SimulationBoundary()
   {
   }

   void SimulationBoundary::convert(const std::map<std::string,int>& bcIds)
   {
      for(auto mapIt = bcIds.cbegin(); mapIt != bcIds.cend(); ++mapIt)
      {
         std::size_t pn = Hasher::makeId(mapIt->first);

         this->mBcs.insert(std::make_pair(pn, mapIt->second));
      }

   }

   std::map<std::string,int>  SimulationBoundary::getTagMap() const
   {
      std::map<std::string,int>  tagMap;

      for(auto mapIt = this->mBcs.cbegin(); mapIt != this->mBcs.cend(); ++mapIt)
      {
         tagMap.insert(std::make_pair(PhysicalNames::Coordinator::tag(mapIt->first), mapIt->second));
      }

      return tagMap;
   }

   int  SimulationBoundary::bcId(const std::size_t id) const
   {
      assert(this->mBcs.count(id) > 0);
      return this->mBcs.at(id);
   }

}
