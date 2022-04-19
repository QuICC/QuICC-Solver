/** 
 * @file Io.hpp 
 * @brief Implementation of the IO node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_IO_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_IO_HPP

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

namespace Framework {

   /**
    * @brief Implementation of the IO node of the configuration file
    */
   class Io: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         Io();

         /**
          * @brief Destructor
          */
         virtual ~Io();

         /**
          * @brief Check compatibility of data
          *
          * \mhdBug Check is not yet implemented
          */
         virtual void checkData();
         
      protected:
         /**
          * @brief Tag name of the parent node
          */
         static const std::string  PARENTTAG;

         /**
          * @brief Initialise component
          */
         void init();

      private:
   };

   /// Typedef for a shared pointer of a IO node
   typedef std::shared_ptr<Io> SharedIo;

   /// Typedef for a const shared pointer of a IO node
   typedef std::shared_ptr<const Io> SharedCIo;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_IO_HPP
