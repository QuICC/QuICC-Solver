/** 
 * @file Block.cpp
 * @brief Source of the model configuration block
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/Model/Block.hpp"

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/Io/Config/Model/Io.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Model {

   const IConfigurationNode& Block::node(const std::string tag) const
   {
      auto id = Hasher::makeId(tag);

      return IConfigurationBlock::node(static_cast<NodeMap::key_type>(id));
   }

   void Block::addNodes(const std::map<std::string,std::map<std::string,int> >& nodes)
   {
      std::string tag;
      NodeMap::key_type id;
      std::shared_ptr<IConfigurationNode> spN;

      for(auto&& [tag, opts]: nodes)
      {
         id = Hasher::makeId(tag);
         spN = std::make_shared<Io>(tag, opts);
         IConfigurationBlock::addNode(id, spN);
      }
   }

   const std::string Block::TAG = "model";

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
