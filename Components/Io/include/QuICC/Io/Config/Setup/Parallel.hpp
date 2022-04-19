/**
 * @file Parallel.hpp
 * @brief Implementation of the setup parallel node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_SETUP_PARALLEL_HPP
#define QUICC_IO_CONFIG_SETUP_PARALLEL_HPP

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

namespace Setup {

   /**
    * @brief Implementation of the setup parallel node of the configuration file
    */
   class Parallel: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit Parallel();

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
          *
          * @param names Names of the parameters
          */
         void init();

      private:
   };

   /// Typedef for a shared pointer of a setup parallel node
   typedef std::shared_ptr<Parallel> SharedParallel;

   /// Typedef for a const shared pointer of a setup parallel node
   typedef std::shared_ptr<const Parallel> SharedCParallel;

}
}
}
}

#endif // QUICC_IO_CONFIG_SETUP_PARALLEL_HPP
