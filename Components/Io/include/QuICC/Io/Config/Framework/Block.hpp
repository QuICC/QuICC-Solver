/**
 * @file Block.hpp
 * @brief Implementation of the setup configuration block
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_BLOCK_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_BLOCK_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Io/Config/IConfigurationNode.hpp"
#include "QuICC/Io/Config/IConfigurationBlock.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Framework {

   /**
    * @name Enum for framework configuration block
    */
   enum NodeId {TRUNCATION, PARALLEL, TIMESTEPPING, RUN, IO, STATISTICS};

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
         void init(const int dim, const std::vector<bool>& isPeriodicBox);

         /**
          * @brief Get boundary node
          */
         SharedCIConfigurationNode spNode(NodeId id) const;

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

#endif // QUICC_IO_CONFIG_FRAMEWORK_BLOCK_HPP
