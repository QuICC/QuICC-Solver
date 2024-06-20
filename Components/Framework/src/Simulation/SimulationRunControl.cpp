/**
 * @file SimulationRunControl.cpp
 * @brief Implementation of a general simulation control structure
 */

// System includes
//
#include <csignal>

// External includes
//

// Class include
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Simulation/SimulationRunControl.hpp"

// Project includes
//
#include "QuICC/RuntimeStatus/GoOn.hpp"
#include "QuICC/RuntimeStatus/Stop.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   std::size_t SimulationRunControl::SIGNAL_STATUS = RuntimeStatus::GoOn::id();

   SimulationRunControl::SimulationRunControl()
      : mStatus(RuntimeStatus::GoOn::id()), mCtrlFile(QuICCEnv().allowsIO()), mSteps(0), mMaxSimTime(0.0), mMaxWallTime(0.0)
   {
      // Initialise signal handling
      this->initSignalHandler();
   }

   SimulationRunControl::~SimulationRunControl()
   {
   }

   std::size_t SimulationRunControl::status() const
   {
      return this->mStatus;
   }

   void SimulationRunControl::updateSimulation(const MHDFloat simTime, const MHDFloat simDt)
   {
      // Increment timestep counter
      this->mSteps++;

      // Check for maximum simulation time
      if(this->mMaxSimTime > 0 && simTime > this->mMaxSimTime)
      {
         this->mStatus = RuntimeStatus::Stop::id();

         // Produce a nice looking output to std output
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printCentered(std::cout, "Simulation time limit reached!", '#');
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printNewline(std::cout);
      }

      // Check for maximum simulation steps
      if(this->mMaxSimTime < 0 && this->mSteps >= std::abs(this->mMaxSimTime))
      {
         this->mStatus = RuntimeStatus::Stop::id();

         // Produce a nice looking output to std output
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printCentered(std::cout, "Simulation steps limit reached!", '#');
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printNewline(std::cout);
      }

      // Check if timestepper requested abort due to too small timestep
      if(simDt < 0)
      {
         this->mStatus = RuntimeStatus::Stop::id();

         // Produce a nice looking output to std output
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printCentered(std::cout, "Adaptive timestep failed!", '#');
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printNewline(std::cout);
      }
   }

   void SimulationRunControl::updateCluster(const MHDFloat wallTime)
   {
      // Not ideal but OK for the moment
      this->mCtrlFile.read();

      // Convert wallTime to hours
      MHDFloat hours = wallTime/3600.;

      // Control status
      if(this->mCtrlFile.status() == RuntimeStatus::Stop::id())
      {
         this->mStatus = RuntimeStatus::Stop::id();
      }

      // Check for maximum wall time
      if(this->mMaxWallTime > 0 && hours > this->mMaxWallTime)
      {
         this->mStatus = RuntimeStatus::Stop::id();

         // Produce a nice looking output to std output
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printCentered(std::cout, "Simulation wall time reached!", '#');
         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printNewline(std::cout);
      }

      // Signal status
      if(SimulationRunControl::SIGNAL_STATUS == RuntimeStatus::Stop::id())
      {
         this->mStatus = RuntimeStatus::Stop::id();
      }
   }

   void SimulationRunControl::checkFile()
   {
      // Read input from control file
      this->mCtrlFile.read();
   }

   void SimulationRunControl::setMaxSimTime(const MHDFloat maxTime)
   {
      this->mMaxSimTime = maxTime;
   }

   void SimulationRunControl::setMaxWallTime(const MHDFloat maxTime)
   {
      this->mMaxWallTime = maxTime;
   }

   void SimulationRunControl::printInfo(std::ostream& stream)
   {
      // Create nice looking ouput header
      Tools::Formatter::printNewline(stream);
      Tools::Formatter::printLine(stream, '-');
      Tools::Formatter::printCentered(stream, "Simulation run information", '*');
      Tools::Formatter::printLine(stream, '-');

      std::stringstream oss;

      // get a nice base for info
      int base = 20;
      oss << std::fixed << std::setprecision(1) << this->mSteps;
      base += oss.str().size() + 1;
      oss.str("");

      // Output number of timesteps
      oss << "Timesteps: " << std::fixed << std::setprecision(1) << this->mSteps;
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      Tools::Formatter::printLine(stream, '*');
   }

   void SimulationRunControl::handleSignal(int signum)
   {
      if(signum == SIGUSR1)
      {
         SimulationRunControl::SIGNAL_STATUS = RuntimeStatus::Stop::id();

         Tools::Formatter::printLine(std::cout, '#');
         Tools::Formatter::printCentered(std::cout, "Simulation received stop signal!", '#');
         Tools::Formatter::printLine(std::cout, '#');
      }
   }

   void SimulationRunControl::initSignalHandler()
   {
      signal(SIGUSR1, SimulationRunControl::handleSignal);
   }

}
