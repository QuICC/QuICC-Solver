/** 
 * @file Block.cpp
 * @brief Source of the setup configuration block
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Framework/Block.hpp"

// Project includes
//
#include "QuICC/Io/Config/Framework/Io.hpp"
#include "QuICC/Io/Config/Framework/Parallel.hpp"
#include "QuICC/Io/Config/Framework/Run.hpp"
#include "QuICC/Io/Config/Framework/Statistics.hpp"
#include "QuICC/Io/Config/Framework/Timestepping.hpp"
#include "QuICC/Io/Config/Framework/Truncation.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   SharedCIConfigurationNode Block::spNode(NodeId id) const
   {
      return IConfigurationBlock::spNode(static_cast<NodeMap::key_type>(id));
   }

   const std::string Block::TAG = "framework";

   Block::Block()
      : IConfigurationBlock(Block::TAG)
   {
   }

   Block::~Block()
   {
   }

   void Block::init(const int dim, const std::vector<bool>& isPeriodicBox)
   {
      // Create shared pointer
      SharedIConfigurationNode spNode;

      //
      // Create setup content

      // Add node
      spNode = std::make_shared<Io>();
      this->addNode(static_cast<NodeMap::key_type>(IO), spNode);

      // Add node
      spNode = std::make_shared<Parallel>();
      this->addNode(static_cast<NodeMap::key_type>(PARALLEL), spNode);

      // Add node
      spNode = std::make_shared<Run>();
      this->addNode(static_cast<NodeMap::key_type>(RUN), spNode);

      // Add node
      spNode = std::make_shared<Statistics>();
      this->addNode(static_cast<NodeMap::key_type>(STATISTICS), spNode);

      // Add node
      spNode = std::make_shared<Timestepping>();
      this->addNode(static_cast<NodeMap::key_type>(TIMESTEPPING), spNode);

      // Add node
      spNode = std::make_shared<Truncation>(dim, isPeriodicBox);
      this->addNode(static_cast<NodeMap::key_type>(TRUNCATION), spNode);
   }
}
}
}
}
