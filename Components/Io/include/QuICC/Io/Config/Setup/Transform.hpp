/**
 * @file Transform.hpp
 * @brief Implementation of the setup transform node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_SETUP_TRANSFORM_HPP
#define QUICC_IO_CONFIG_SETUP_TRANSFORM_HPP

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
    * @brief Implementation of the setup transform node of the configuration file
    */
   class Transform: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit Transform(const int dim);

         /**
          * @brief Destructor
          */
         virtual ~Transform();

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
         void init(const int dim);

      private:
   };

   /// Typedef for a shared pointer of a setup transform node
   typedef std::shared_ptr<Transform> SharedTransform;

   /// Typedef for a const shared pointer of a setup transform node
   typedef std::shared_ptr<const Transform> SharedCTransform;

}
}
}
}

#endif // QUICC_IO_CONFIG_SETUP_TRANSFORM_HPP
