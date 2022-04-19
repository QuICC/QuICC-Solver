/**
 * @file Io.hpp
 * @brief Implementation of the a simple IO node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_MODEL_IO_HPP
#define QUICC_IO_CONFIG_MODEL_IO_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Config/IConfigurationNode.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Model {

   /**
    * @brief Implementation of the setup boundary node of the configuration file
    */
   class Io: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit Io(const std::string tag, const std::map<std::string,int>& options);

         /**
          * @brief Destructor
          */
         virtual ~Io();

         /**
          * @brief Check compatibility of data
          */
         virtual void checkData();

      protected:
         /**
          * @brief Initialise component
          */
         void init(const std::map<std::string,int>& options);

      private:
   };

}
}
}
}

#endif // QUICC_IO_CONFIG_MODEL_IO_HPP
