/** 
 * @file Statistics.hpp 
 * @brief Implementation of the statistics node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_STATISTICS_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_STATISTICS_HPP

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
    * @brief Implementation of the statistics node of the configuration file
    */
   class Statistics: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         Statistics();

         /**
          * @brief Destructor
          */
         virtual ~Statistics();

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
   typedef std::shared_ptr<Statistics> SharedStatistics;

   /// Typedef for a const shared pointer of a IO node
   typedef std::shared_ptr<const Statistics> SharedCStatistics;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_STATISTICS_HPP
