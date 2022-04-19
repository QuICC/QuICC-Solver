/**
 * @file Common.hpp
 * @brief Definition of Common component version information
 */

#ifndef QUICC_VERSION_COMMON_HPP
#define QUICC_VERSION_COMMON_HPP

#define QUICC_VERSION_COMMON_MAJOR 1
#define QUICC_VERSION_COMMON_MINOR 0
#define QUICC_VERSION_COMMON_PATCH 0

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
   class Common 
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

#endif // QUICC_VERSION_COMMON_HPP
