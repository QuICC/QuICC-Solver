/** 
 * @file Timestepping.hpp 
 * @brief Implementation of the timestepping node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_TIMESTEPPING_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_TIMESTEPPING_HPP

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
    * @brief Implementation of the timestepping node of the configuration file
    */
   class Timestepping: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         Timestepping();

         /**
          * @brief Destructor
          */
         virtual ~Timestepping();

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

   /// Typedef for a shared pointer of a timestepping node
   typedef std::shared_ptr<Timestepping> SharedTimestepping;

   /// Typedef for a const shared pointer of a timestepping node
   typedef std::shared_ptr<const Timestepping> SharedCTimestepping;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_TIMESTEPPING_HPP
