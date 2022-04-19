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
#include "QuICC/Io/Config/Simulation/Block.hpp"

// Project includes
//
#include "QuICC/Io/Config/Simulation/Boundary.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Simulation {

   SharedCIConfigurationNode Block::spNode(NodeId id) const
   {
      return IConfigurationBlock::spNode(static_cast<NodeMap::key_type>(id));
   }

   void Block::addNode(NodeId id, SharedIConfigurationNode spNode)
   {
      IConfigurationBlock::addNode(static_cast<NodeMap::key_type>(id), spNode);
   }

   const std::string Block::TAG = "simulation";

   Block::Block()
      : IConfigurationBlock(Block::TAG)
   {
   }

   Block::~Block()
   {
   }

   void Block::init()
   {
   }
}
}
}
}
