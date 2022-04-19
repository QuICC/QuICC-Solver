/**
 * @file Framework.hpp
 * @brief Definition of Framework component version information
 */

#ifndef QUICC_VERSION_FRAMEWORK_HPP
#define QUICC_VERSION_FRAMEWORK_HPP

#define QUICC_VERSION_FRAMEWORK_MAJOR 1
#define QUICC_VERSION_FRAMEWORK_MINOR 1
#define QUICC_VERSION_FRAMEWORK_PATCH 0

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Version {

   /**
    * @brief Small class storing version information
    */
   class Framework 
   {
      public:
         /**
          * @brief Major version number, ie MAJOR.MINOR.PATCH
          */
         static const int MAJOR;

         /**
          * @brief Minor version number, ie MAJOR.MINOR.PATCH
          */
         static const int MINOR;

         /**
          * @brief Patch version number, ie MAJOR.MINOR.PATCH
          */
         static const int PATCH;

         /**
          * @brief Get version number string, ie MAJOR.MINOR.PATCH
          */
         static std::string version();
   };

}
}

#endif // QUICC_VERSION_FRAMEWORK_HPP
