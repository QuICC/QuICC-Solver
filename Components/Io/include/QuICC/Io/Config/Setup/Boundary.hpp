/**
 * @file Boundary.hpp
 * @brief Implementation of the setup boundary node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_SETUP_BOUNDARY_HPP
#define QUICC_IO_CONFIG_SETUP_BOUNDARY_HPP

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
    * @brief Implementation of the setup boundary node of the configuration file
    */
   class Boundary: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit Boundary();

         /**
          * @brief Destructor
          */
         virtual ~Boundary();

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

   /// Typedef for a shared pointer of a boundary node
   typedef std::shared_ptr<Boundary> SharedBoundary;

   /// Typedef for a const shared pointer of a boundary node
   typedef std::shared_ptr<const Boundary> SharedCBoundary;

}
}
}
}

#endif // QUICC_IO_CONFIG_SETUP_BOUNDARY_HPP
