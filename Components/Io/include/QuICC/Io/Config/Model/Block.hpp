/**
 * @file Block.hpp
 * @brief Implementation of the model configuration block
 */

#ifndef QUICC_IO_CONFIG_MODEL_BLOCK_HPP
#define QUICC_IO_CONFIG_MODEL_BLOCK_HPP

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

namespace Model {

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
          * @brief Initialise the different parts of the model part
          */
         void init();

         /**
          * @brief Get boundary node
          */
         const IConfigurationNode& node(const std::string tag) const;

         /**
          * @brief Add model configuration
          */
         void addNodes(const std::map<std::string,std::map<std::string,int> >& nodes);

      protected:

      private:
         /**
          * @brief Name of the model XML tag
          */
         static const std::string TAG;
   };

}
}
}
}

#endif // QUICC_IO_CONFIG_MODEL_BLOCK_HPP
