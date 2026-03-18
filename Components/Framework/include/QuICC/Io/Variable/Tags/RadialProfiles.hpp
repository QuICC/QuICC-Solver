/** 
 * @file RadialProfiles.hpp
 * @brief Definitions and names use by the Radial Profiles writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_RADIALPROFILES
#define QUICC_IO_VARIABLE_TAGS_RADIALPROFILES

// System includes
//
#include <string>

// Project includes
//

namespace QuICC {

namespace Io {

namespace Variable {

namespace Tags {

   /**
    * @brief Definitions and names use by the RadialProfiles writer
    */
   class RadialProfiles
   {
      public:
         /**
          * @brief HEADER part for RadialProfiles file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for RadialProfiles file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of RadialProfiles file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of RadialProfiles file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         RadialProfiles();

         /**
         * @brief Destructor
         */
         ~RadialProfiles();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_RADIALPROFILES
