/**
 * @file SimulationConfig.cpp
 * @brief Source of the implementation of the IO control system
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Simulation/SimulationConfig.hpp"

// Project includes
//
#include "QuICC/Enums/SplittingTools.hpp"
#include "QuICC/Timestep/Id/Coordinator.hpp"
#include "QuICC/Timestep/Id/registerAll.hpp"

namespace QuICC {

   SimulationConfig::SimulationConfig(Io::Config::SharedConfigurationReader spCfgFile)
      : mspCfgFile(spCfgFile)
   {
   }

   SimulationConfig::~SimulationConfig()
   {
   }

   ArrayI SimulationConfig::dimension() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      // Get truncation map
      std::map<std::string,int>  trunc = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::TRUNCATION)->iTags().map();

      // Create storage for the dimensions
      ArrayI dim(trunc.size());

      // Extrac dimension from truncation read from file
      int i = 0;
      for(auto itI = trunc.cbegin(); itI != trunc.cend(); itI++)
      {
         dim(i) = itI->second;
         i++;
      }

      return dim;
   }

   ArrayI SimulationConfig::transformSetup() const
   {
      auto nodeId = Io::Config::Setup::TRANSFORM;

      // Get map
      int dim = this->mspCfgFile->spSetup()->spNode(nodeId)->iTags().size();

      ArrayI impl(dim);
      int i = 0;
      auto range = this->mspCfgFile->spSetup()->spNode(nodeId)->iTags().crange();
      for(auto it = range.first; it != range.second; ++it)
      {
         impl(i) = it->second;
         i++;
      }

      return impl;
   }

   ArrayI SimulationConfig::boundarySetup() const
   {
      auto nodeId = Io::Config::Setup::BOUNDARY;

      // Get map
      int dim = this->mspCfgFile->spSetup()->spNode(nodeId)->iTags().size();

      ArrayI tags(dim);
      int i = 0;
      auto range = this->mspCfgFile->spSetup()->spNode(nodeId)->iTags().crange();
      for(auto it = range.first; it != range.second; ++it)
      {
         tags(i) = it->second;
         i++;
      }

      return tags;
   }

   Array SimulationConfig::boxScale() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto cfgTrunc = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::TRUNCATION);
      // Create storage for box scales
      Array box = Array::Ones(cfgTrunc->iTags().size());

      // Get box scale truncation map
      std::map<std::string,MHDFloat>  trunc = cfgTrunc->fTags().map();

      // Extract box scale from truncation read from file
      for(auto it = trunc.cbegin(); it != trunc.cend(); it++)
      {
         if(it->first == "kc1D")
         {
            box(0) = it->second/trunc.find("box1D")->second;
         } else if(it->first == "kc2D")
         {
            box(1) = it->second/trunc.find("box2D")->second;
         } else if(it->first == "kc3D")
         {
            box(2) = it->second/trunc.find("box3D")->second;
         }
      }

      return box;
   }

   int SimulationConfig::nCpu() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::PARALLEL)->iTags().value("cpus");
   }

   Splitting::Algorithms::Id SimulationConfig::algorithm() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto tag = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::PARALLEL)->sTags().value("algorithm");
      auto id = Splitting::getAlgorithmId(tag);

      return id;
   }

   Splitting::Groupers::Id SimulationConfig::grouper() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto tag = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::PARALLEL)->sTags().value("grouper");
      auto id = Splitting::getGrouperId(tag);

      return id;
   }

   const std::map<std::string, MHDFloat>& SimulationConfig::physical() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spSimulation()->spNode(Io::Config::Simulation::PHYSICAL)->fTags().map();
   }

   std::map<std::string, MHDFloat>& SimulationConfig::rPhysical()
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->rspSimulation()->rspNode(Io::Config::Simulation::PHYSICAL)->fTags().rMap();
   }

   const std::map<std::string, int>& SimulationConfig::boundary() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->spSimulation()->spNode(Io::Config::Simulation::BOUNDARY)->iTags().map();
   }

   const std::map<std::string, int>& SimulationConfig::model(const std::string tag) const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      return this->mspCfgFile->model().node(tag).iTags().map();
   }

   Array SimulationConfig::run() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      Array run(2);

      auto spRun = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::RUN);
      // Get simulation time configuration
      run(0) = spRun->fTags().value("sim");

      // Get wall time configuration
      run(1) = spRun->fTags().value("wall");

      return run;
   }

   Array SimulationConfig::timestepping() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      Array tstep(3);

      // Get timestepping time configuration
      auto spTime = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::TIMESTEPPING);
      tstep(0) = spTime->fTags().value("time");

      // Get timestepping timestep configuration
      tstep(1) = spTime->fTags().value("timestep");

      // Get timestepping error configuration
      tstep(2) = spTime->fTags().value("error");

      return tstep;
   }

   std::size_t SimulationConfig::timestepper() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto s = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::TIMESTEPPING)->sTags().value("scheme");
      std::string tag = s;
      std::transform(s.cbegin(), s.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

      // Register all timestepper IDs
      Timestep::Id::registerAll();

      std::size_t id = 0;
      for(auto&& e: Timestep::Id::Coordinator::map())
      {
         if(Timestep::Id::Coordinator::tag(e.first) == tag)
         {
            id = e.first;
            break;
         }
      }

      if(id == 0)
      {
         throw std::logic_error("timestepper scheme was not recognized");
      }

      return id;
   }
}
