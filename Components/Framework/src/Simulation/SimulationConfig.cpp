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
#include "QuICC/Timestep/Id/registerAllPc.hpp"
#include "QuICC/Timestep/Id/registerAllRkcb.hpp"
#include "QuICC/Transform/Setup/registerAll.hpp"
#include "QuICC/Bc/Scheme/Coordinator.hpp"
#include "QuICC/Bc/Scheme/registerAll.hpp"
#include "QuICC/Bc/Name/Coordinator.hpp"
#include "QuICC/Bc/Name/registerAll.hpp"
#include "QuICC/Tag/Generic/registerAll.hpp"

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

   std::map<std::size_t, std::vector<std::size_t>> SimulationConfig::transformSetup() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      std::map<std::size_t,std::vector<std::size_t>> impl;

      const auto nodeId = Io::Config::Setup::TRANSFORM;

      // Register all transform implementation IDs
      Transform::Setup::registerAll();

      int i = 0;
      auto range = this->mspCfgFile->spSetup()->spNode(nodeId)->sTags().crange();
      for(auto sit = range.first; sit != range.second; ++sit)
      {
         auto it = impl.emplace(i, std::vector<std::size_t>());

         std::string tag = sit->second;
         std::transform(sit->second.cbegin(), sit->second.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

         if(tag.size() > 0)
         {
            std::string tag_;
            std::size_t start = 0;
            // Extract individual tags from comma separated list
            do
            {
               std::size_t pos = tag.find(',', start);
               if(pos == std::string::npos)
               {
                  pos = tag.size();
               }
               tag_ = tag.substr(start, pos-start);
               start = pos + 1;

               // Get ID for tag
               std::size_t id = 0;
               for(auto&& e: Transform::Setup::Coordinator::map())
               {
                  if(Transform::Setup::Coordinator::tag(e.first) == tag_)
                  {
                     id = e.first;
                     break;
                  }
               }

               if(id == 0)
               {
                  throw std::logic_error("transform setup was not recognized");
               }

               it.first->second.push_back(id);
            }
            while(start < tag.size()-1);
         }
         i++;
      }

      return impl;
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

   std::list<int>  SimulationConfig::cpuFactors() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto tag = this->mspCfgFile->spFramework()->spNode(Io::Config::Framework::PARALLEL)->sTags().value("decomposition");
      std::transform(tag.cbegin(), tag.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

      std::list<int> factorList;
      if(tag != "auto" && tag.size() > 0)
      {
         std::string tag_;
         std::size_t start = 0;
         // Extract individual factors from comma separated list
         do
         {
            std::size_t pos = tag.find(',', start);
            if(pos == std::string::npos)
            {
               pos = tag.size();
            }
            tag_ = tag.substr(start, pos-start);
            start = pos + 1;

            // Get ID for tag
            int f;
            try
            {
               f = std::stoi(tag_);
            }
            catch(const std::invalid_argument& e)
            {
               throw std::logic_error("Decomposition format is wrong. It should be 'auto' or a comma separated list of factors");
            }
            factorList.emplace_back(f);
         }
         while(start < tag.size());
      }

      return factorList;
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

   std::map<std::string, std::size_t> SimulationConfig::boundary() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto s = this->mspCfgFile->spSimulation()->spNode(Io::Config::Simulation::BOUNDARY)->sTags().map();
      std::map<std::string, std::size_t> bcMap;

      // Register all boundary condition names
      Bc::Name::registerAll();

      for(auto bc: s)
      {
         auto s = bc.second;
         std::string tag = s;
         std::transform(s.cbegin(), s.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

         std::size_t id = 0;
         for(auto&& e: Bc::Name::Coordinator::map())
         {
            if(Bc::Name::Coordinator::tag(e.first) == tag)
            {
               id = e.first;
               break;
            }
         }

         if(id == 0)
         {
            throw std::logic_error("boundary condition was not recognized");
         }

         bcMap.emplace(std::make_pair(bc.first, id));
      }

      return bcMap;
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
      Timestep::Id::registerAllPc();
      Timestep::Id::registerAllRkcb();

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

      // For the time being only the Predictor-correct scheme supports split equations
      if(this->splitEquation() && (id != Timestep::Id::ImexPc2::id()))
      {
         throw std::logic_error("Split equations are only supported with PC2 scheme");
      }

      return id;
   }

   std::size_t SimulationConfig::bcScheme() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto s = this->mspCfgFile->spSetup()->spNode(Io::Config::Setup::BOUNDARY)->sTags().value("scheme");
      std::string tag = s;
      std::transform(s.cbegin(), s.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

      // Register all boundary condition scheme IDs
      Bc::Scheme::registerAll();

      std::size_t id = 0;
      for(auto&& e: Bc::Scheme::Coordinator::map())
      {
         if(Bc::Scheme::Coordinator::tag(e.first) == tag)
         {
            id = e.first;
            break;
         }
      }

      if(id == 0)
      {
         throw std::logic_error("boundary condition scheme was not recognized");
      }

      return id;
   }

   bool SimulationConfig::splitEquation() const
   {
      // Safety assert for non NULL pointer
      assert(this->mspCfgFile);

      auto s = this->mspCfgFile->spSetup()->spNode(Io::Config::Setup::MODEL)->sTags().value("split_equation");
      std::string tag = s;
      std::transform(s.cbegin(), s.cend(), tag.begin(), [](unsigned char c) { return std::tolower(c); });

      // Register all boundary condition scheme IDs
      Tag::Generic::registerAll();

      std::size_t id = 0;
      for(auto&& e: Tag::Generic::Coordinator::map())
      {
         if(Tag::Generic::Coordinator::tag(e.first) == tag)
         {
            id = e.first;
            break;
         }
      }

      if(id == 0)
      {
         throw std::logic_error("on/off tag for split_equation was not recognized");
      }

      bool enable = (Tag::Generic::On::id() == id);

      return enable;
   }
}
