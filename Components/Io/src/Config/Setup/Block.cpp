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
#include "QuICC/Io/Config/Setup/Block.hpp"

// Project includes
//
#include "QuICC/Io/Config/Setup/Boundary.hpp"
#include "QuICC/Io/Config/Setup/Model.hpp"
#include "QuICC/Io/Config/Setup/Parallel.hpp"
#include "QuICC/Io/Config/Setup/Timestepper.hpp"
#include "QuICC/Io/Config/Setup/Transform.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Setup {

   SharedCIConfigurationNode Block::spNode(NodeId id) const
   {
      return IConfigurationBlock::spNode(static_cast<NodeMap::key_type>(id));
   }

   const std::string Block::TAG = "setup";

   Block::Block()
      : IConfigurationBlock(Block::TAG)
   {
   }

   Block::~Block()
   {
   }

   void Block::init(const int dim)
   {
      // Create shared pointer
      SharedIConfigurationNode spNode;

      //
      // Create setup content

      // Add model part
      spNode = std::make_shared<Model>();
      this->addNode(static_cast<NodeMap::key_type>(MODEL), spNode);

      // Add boundary part
      spNode = std::make_shared<Boundary>();
      this->addNode(static_cast<NodeMap::key_type>(BOUNDARY), spNode);

      // Add parallel part
      spNode = std::make_shared<Parallel>();
      this->addNode(static_cast<NodeMap::key_type>(PARALLEL), spNode);

      // Add timestepper part
      spNode = std::make_shared<Timestepper>();
      this->addNode(static_cast<NodeMap::key_type>(TIMESTEPPER), spNode);

      // Add transform part
      spNode = std::make_shared<Transform>(dim);
      this->addNode(static_cast<NodeMap::key_type>(TRANSFORM), spNode);
   }
}
}
}
}
