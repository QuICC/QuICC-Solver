/**
 * @file Model.hpp
 * @brief Implementation of the setup model node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_SETUP_MODEL_HPP
#define QUICC_IO_CONFIG_SETUP_MODEL_HPP

// System includes
//
#include <string>
#include <memory>

// Project includes
//
#include "QuICC/Io/Config/IConfigurationNode.hpp"

namespace QuICC {

namespace Io {

namespace Config {

namespace Setup {

   /**
    * @brief Implementation of the setup model node of the configuration file
    */
   class Model: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit Model();

         /**
          * @brief Destructor
          */
         virtual ~Model() = default;

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

   /// Typedef for a shared pointer of a model node
   typedef std::shared_ptr<Model> SharedModel;

   /// Typedef for a const shared pointer of a model node
   typedef std::shared_ptr<const Model> SharedCModel;

}
}
}
}

#endif // QUICC_IO_CONFIG_SETUP_MODEL_HPP
