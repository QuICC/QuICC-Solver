/**
 * @file SimulationIoTools.cpp
 * @brief Source of the implementation of the tools for IO related calculations
 */

// System includes
//

// Project includes
//
#include "QuICC/Simulation/SimulationIoTools.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

   SimulationIoTools::SimulationIoTools()
   {
   }

   SimulationIoTools::~SimulationIoTools()
   {
   }

   void SimulationIoTools::updateHeavyAscii(SimulationIoTools::ascii_iterator asciiBegin,  SimulationIoTools::ascii_iterator asciiEnd, Transform::TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("SimulationIoTools::updateHeavyAscii");

      ascii_iterator it;
      for(it = asciiBegin; it != asciiEnd; ++it)
      {
         if((*it)->isHeavy() && (*it)->isActive())
         {
            (*it)->compute(coord);
         }
      }
   }

   void SimulationIoTools::updateStatsPre(SimulationIoTools::stats_iterator statsBegin,  SimulationIoTools::stats_iterator statsEnd, Transform::TransformCoordinatorType& coord)
   {
      stats_iterator it;
      for(it = statsBegin; it != statsEnd; ++it)
      {
         (*it)->preCompute(coord);
      }
   }

   void SimulationIoTools::updateStats(SimulationIoTools::stats_iterator statsBegin,  SimulationIoTools::stats_iterator statsEnd, Transform::TransformCoordinatorType& coord)
   {
      stats_iterator it;
      for(it = statsBegin; it != statsEnd; ++it)
      {
         (*it)->compute(coord);
      }
   }

   void SimulationIoTools::updateStatsPost(SimulationIoTools::stats_iterator statsBegin,  SimulationIoTools::stats_iterator statsEnd, Transform::TransformCoordinatorType& coord)
   {
      stats_iterator it;
      for(it = statsBegin; it != statsEnd; ++it)
      {
         (*it)->postCompute(coord);
      }
   }
}
