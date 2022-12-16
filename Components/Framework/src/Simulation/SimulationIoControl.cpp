/**
 * @file SimulationIoControl.cpp
 * @brief Source of the implementation of the IO control system
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Simulation/SimulationIoControl.hpp"

// Project includes
//
#include "QuICC/Version/Common.hpp"
#include "QuICC/Version/Framework.hpp"
#include "QuICC/Version/Io.hpp"
#include "QuICC/Version/Polynomial.hpp"
#include "QuICC/Version/PyQuICC.hpp"
#include "QuICC/Version/SparseSM.hpp"
#include "QuICC/Version/Transform.hpp"

namespace QuICC {

   SimulationIoControl::SimulationIoControl()
      : mSteps(0), mAsciiRate(-1), mHdf5Rate(-1), mStatsRate(-1), mStatsAvgRate(-1), mActiveStatsUpdate(false), mActiveStatsWrite(true)
   {
   }

   SimulationIoControl::~SimulationIoControl()
   {
      // Finalize ASCII and HDF5 writers
      this->finalizeWriters();
   }

   void SimulationIoControl::setConfigurationFile(Io::Config::SharedConfigurationReader spCfgFile)
   {
      this->mspCfgFile = spCfgFile;

      this->mspCfg = std::make_shared<SimulationConfig>(this->mspCfgFile);
   }

   void SimulationIoControl::init()
   {
      // Init configuration file
      this->initCfg();

      // Init StdOutPipe output file
      this->initStdOut();

      // Print parameters file parameters
      this->mspCfgFile->printInfo();

      // Print version information
      Tools::Formatter::printNewline(std::cout);
      Tools::Formatter::printLine(std::cout, '=');
      Tools::Formatter::printCentered(std::cout, "QuICC", '*');
      Tools::Formatter::printLine(std::cout, '-');
      Tools::Formatter::printCentered(std::cout, "Common: " + Version::Common::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "Framework: " + Version::Framework::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "Io: " + Version::Io::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "Polynomial: " + Version::Polynomial::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "PyQuICC: " + Version::PyQuICC::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "SparseSM: " + Version::SparseSM::version(), ' ');
      Tools::Formatter::printCentered(std::cout, "Transform: " + Version::Transform::version(), ' ');
      Tools::Formatter::printLine(std::cout, '=');
      Tools::Formatter::printNewline(std::cout);
   }

   void SimulationIoControl::cleanup()
   {
      // reset the configuration file
      this->mspCfgFile.reset();
   }

   void SimulationIoControl::update()
   {
      // Increment timestep counter
      this->mSteps++;
   }

   bool SimulationIoControl::isAsciiTime() const
   {
      return (this->mAsciiRate > 0 && this->mSteps % this->mAsciiRate == 0);
   }

   bool SimulationIoControl::isHdf5Time() const
   {
      return (this->mHdf5Rate > 0 && this->mSteps % this->mHdf5Rate == 0);
   }

   void SimulationIoControl::activateStats()
   {
      bool writeTrigger = (this->mStatsRate > 0 && this->mSteps % this->mStatsRate == 0);
      bool updateTrigger = (this->mStatsAvgRate > 0 && this->mSteps % this->mStatsAvgRate == 0);

      this->mActiveStatsUpdate = (writeTrigger || updateTrigger);
   }

   void SimulationIoControl::disableStats()
   {
      this->mActiveStatsUpdate = false;
   }

   bool SimulationIoControl::isStatsTime() const
   {
      bool writeTrigger = (this->mStatsRate > 0 && this->mSteps % this->mStatsRate == 0);
      return writeTrigger;
   }

   bool SimulationIoControl::isStatsUpdateTime() const
   {
      return this->mActiveStatsUpdate;
   }

   void SimulationIoControl::writeFiles(const MHDFloat time, const MHDFloat timestep)
   {
      if(this->isAsciiTime())
      {
         this->writeAscii(time,timestep);
      }

      if(this->isHdf5Time())
      {
         this->writeHdf5(time,timestep);
      }

      // Activate stats
      this->activateStats();

      if(this->isStatsTime())
      {
         this->prepareStats(time,timestep);
         this->mActiveStatsWrite = true;
      }
   }

   void SimulationIoControl::addAsciiOutputFile(Io::Variable::SharedIVariableAsciiWriter spOutFile)
   {
      this->mAsciiWriters.push_back(spOutFile);
   }

   void SimulationIoControl::addHdf5OutputFile(Io::Variable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mHdf5Writers.push_back(spOutFile);
   }

   void SimulationIoControl::addStatsOutputFile(Io::Stats::SharedIStatisticsAsciiWriter spOutFile)
   {
      this->mStatsWriters.push_back(spOutFile);
   }

   void SimulationIoControl::initWriters()
   {
      // First check that all ASCII writers are full
      SimulationIoControl::ascii_iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         if(!(*itAscii)->isFull())
         {
            throw std::logic_error("There are missing variables in the ASCII writers");
         }
      }

      // First check that all HDF5 writers are full
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         if(!(*itHdf5)->isFull())
         {
            throw std::logic_error("There are missing variables in the HDF5 writers");
         }
      }

      // First check that all ASCII writers are full
      SimulationIoControl::stats_iterator itStats;
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         if(!(*itStats)->isFull())
         {
            throw std::logic_error("There are missing variables in the statistics writers");
         }
      }

      // Create physical parameters map
      std::map<std::string,MHDFloat> phys = this->config().physical();
      std::map<std::string,std::size_t> boundary = this->config().boundary();

      // Iterate over all ASCII writer
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->setPhysical(phys, boundary);
         (*itAscii)->init();
      }

      // Iterate over all HDF5 writer
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->setPhysical(phys, boundary);
         (*itHdf5)->init();
      }

      // Iterate over all statistics writer
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         (*itStats)->setPhysical(phys, boundary);
         (*itStats)->init();
      }
   }

   void SimulationIoControl::finalizeWriters()
   {
      // Iterate over all ASCII writer
      SimulationIoControl::ascii_iterator itAscii;
      for(itAscii = this->mAsciiWriters.begin(); itAscii < this->mAsciiWriters.end(); itAscii++)
      {
         (*itAscii)->finalize();
      }
      this->mAsciiWriters.clear();

      // Iterate over all HDF5 writer
      SimulationIoControl::hdf5_iterator itHdf5;
      for(itHdf5 = this->mHdf5Writers.begin(); itHdf5 < this->mHdf5Writers.end(); itHdf5++)
      {
         (*itHdf5)->finalize();
      }
      this->mHdf5Writers.clear();

      // Iterate over all statistics writer
      SimulationIoControl::stats_iterator itStats;
      for(itStats = this->mStatsWriters.begin(); itStats < this->mStatsWriters.end(); itStats++)
      {
         (*itStats)->finalize();
      }
      this->mStatsWriters.clear();
   }

   void SimulationIoControl::writeAscii(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all ASCII writer
      SimulationIoControl::ascii_iterator it;
      for(it = this->mAsciiWriters.begin(); it < this->mAsciiWriters.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
         (*it)->write();
      }
   }

   void SimulationIoControl::writeHdf5(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all HDF5 writer
      SimulationIoControl::hdf5_iterator it;
      for(it = this->mHdf5Writers.begin(); it < this->mHdf5Writers.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
         (*it)->write();
      }
   }

   void SimulationIoControl::prepareStats(const MHDFloat time, const MHDFloat timestep)
   {
      // Iterate over all statistics writer
      SimulationIoControl::stats_iterator it;
      for(it = this->mStatsWriters.begin(); it < this->mStatsWriters.end(); it++)
      {
         (*it)->setSimTime(time,timestep);
      }
   }

   void SimulationIoControl::writeStats()
   {
      if(this->mActiveStatsWrite)
      {
         // Iterate over all statistics writer
         SimulationIoControl::stats_iterator it;
         for(it = this->mStatsWriters.begin(); it < this->mStatsWriters.end(); it++)
         {
            (*it)->write();
         }

         this->mActiveStatsWrite = false;
      }
   }

   void SimulationIoControl::initCfg()
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      // Initialise config file
      this->mspCfgFile->init();

      // Read configuration data from config file
      this->mspCfgFile->read();

      // Close the config file
      this->mspCfgFile->finalize();

      // Set ASCII rate from config file
      auto spIo = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::IO);
      this->mAsciiRate = spIo->fTags().value("ascii");

      // Set HDF5 rate from config file
      this->mHdf5Rate = spIo->fTags().value("hdf5");

      // Set statistics rate from config file
      auto spStats = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::STATISTICS);
      this->mStatsRate = spStats->fTags().value("output_rate");
      this->mStatsAvgRate = spStats->fTags().value("time_avg_rate");
   }

   void SimulationIoControl::initStdOut()
   {
      // Create the StdOutPipe ouput
      auto spStd = std::make_shared<Io::Ascii::StdOutPipe>("OUT");

      // Store the shared pointer
      this->mspStdOut = spStd;
   }

   const SimulationConfig& SimulationIoControl::config() const
   {
     return *this->mspCfg;
   }

   SimulationConfig& SimulationIoControl::rConfig()
   {
     return *this->mspCfg;
   }

   SimulationIoControl::ascii_iterator  SimulationIoControl::beginAscii()
   {
      return this->mAsciiWriters.begin();
   }

   SimulationIoControl::ascii_iterator  SimulationIoControl::endAscii()
   {
      return this->mAsciiWriters.end();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::beginHdf5()
   {
      return this->mHdf5Writers.begin();
   }

   SimulationIoControl::hdf5_iterator  SimulationIoControl::endHdf5()
   {
      return this->mHdf5Writers.end();
   }

   SimulationIoControl::stats_iterator  SimulationIoControl::beginStats()
   {
      return this->mStatsWriters.begin();
   }

   SimulationIoControl::stats_iterator  SimulationIoControl::endStats()
   {
      return this->mStatsWriters.end();
   }
}
