/**
 * @file Block.hpp
 * @brief Implementation of the setup configuration block
 */

#ifndef QUICC_IO_CONFIG_SIMULATION_BLOCK_HPP
#define QUICC_IO_CONFIG_SIMULATION_BLOCK_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Io/Config/IConfigurationNode.hpp"
#include "QuICC/Io/Config/IConfigurationBlock.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Simulation {

   /**
    * @name Enum of the possible blocks
    */
   enum NodeId {PHYSICAL, BOUNDARY};

   /**
    * @brief Implementation of the base for a configuration file
    */
   class Block: public IConfigurationBlock
   {
      public:
         /**
          * @brief Constructor
          */
         Block();

         /**
          * @brief Destructor
          */
         virtual ~Block();

         /**
          * @brief Initialise the different parts of the setup part
          */
         void init();

         /**
          * @brief Get boundary node
          */
         SharedCIConfigurationNode spNode(NodeId id) const;

         /**
          * @brief Add configuration part to simulation configuration
          *
          * @param id      ID of the simulation configuration block
          * @param spNode  Configuration node to add
          */
         void addNode(NodeId id, SharedIConfigurationNode spNode);

      protected:

      private:
         /**
          * @brief Name of the setup XML tag
          */
         static const std::string TAG;
   };

}
}
}
}

#endif // QUICC_IO_CONFIG_SIMULATION_BLOCK_HPP
