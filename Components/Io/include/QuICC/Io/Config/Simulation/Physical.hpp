/**
 * @file Physical.hpp
 * @brief Implementation of the physical node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_SIMULATION_PHYSICAL_HPP
#define QUICC_IO_CONFIG_SIMULATION_PHYSICAL_HPP

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

namespace Simulation {

   /**
    * @brief Implementation of the physical node of the configuration file
    */
   class Physical: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          *
          * @param names Names of the parameters
          */
         explicit Physical(const std::vector<std::string>& names);

         /**
          * @brief Destructor
          */
         virtual ~Physical();

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
         void init(const std::vector<std::string>& names);

      private:
   };

   /// Typedef for a shared pointer of a physical node
   typedef std::shared_ptr<Physical> SharedPhysical;

   /// Typedef for a const shared pointer of a physical node
   typedef std::shared_ptr<const Physical> SharedCPhysical;

}
}
}
}

#endif // QUICC_IO_CONFIG_SIMULATION_PHYSICAL_HPP
