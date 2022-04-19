/** 
 * @file Parallel.hpp 
 * @brief Implementation of the parallel node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_PARALLEL_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_PARALLEL_HPP

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
    * @brief Implementation of the parallel node of the configuration file
    */
   class Parallel: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         Parallel();

         /**
          * @brief Destructor
          */
         virtual ~Parallel();

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

   /// Typedef for a shared pointer of a parallel node
   typedef std::shared_ptr<Parallel> SharedParallel;

   /// Typedef for a const shared pointer of a parallel node
   typedef std::shared_ptr<const Parallel> SharedCParallel;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_PARALLEL_HPP
