/** 
 * @file Run.hpp 
 * @brief Implementation of the run node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_RUN_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_RUN_HPP

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
    * @brief Implementation of the run node of the configuration file
    */
   class Run: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         Run();

         /**
          * @brief Destructor
          */
         virtual ~Run();

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

   /// Typedef for a shared pointer of a run node
   typedef std::shared_ptr<Run> SharedRun;

   /// Typedef for a const shared pointer of a run node
   typedef std::shared_ptr<const Run> SharedCRun;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_RUN_HPP
