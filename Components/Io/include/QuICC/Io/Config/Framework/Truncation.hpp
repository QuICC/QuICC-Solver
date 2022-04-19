/** 
 * @file Truncation.hpp 
 * @brief Implementation of the truncation node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_FRAMEWORK_TRUNCATION_HPP
#define QUICC_IO_CONFIG_FRAMEWORK_TRUNCATION_HPP

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
    * @brief Implementation of the truncation node of the configuration file
    *
    * \mhdTodo Strides have been desactivated (need first a clear idea how to use them)
    */
   class Truncation: public IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimensionality of truncation
          */
         explicit Truncation(const int dim, const std::vector<bool>& isPeriodicBox);

         /**
          * @brief Destructor
          */
         virtual ~Truncation();

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
          * @param dim Dimensionality of truncation
          */
         void init(const int dim, const std::vector<bool>& isPeriodicBox);

      private:

   };

   /// Typedef for a shared pointer of a truncation node
   typedef std::shared_ptr<Truncation> SharedTruncation;

   /// Typedef for a const shared pointer of a truncation node
   typedef std::shared_ptr<const Truncation> SharedCTruncation;

}
}
}
}

#endif // QUICC_IO_CONFIG_FRAMEWORK_TRUNCATION_HPP
