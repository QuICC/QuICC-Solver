/** 
 * @file Enstrophy.hpp
 * @brief Definitions and names use by the enstrophy writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_ENSTROPHY
#define QUICC_IO_VARIABLE_TAGS_ENSTROPHY

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Variable {

namespace Tags {

   /**
    * @brief Definitions and names use by the enstrophy writers
    */
   class Enstrophy
   {
      public:
         /**
          * @brief HEADER part for enstrophy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for enstrophy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of enstrophy files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of enstrophy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Enstrophy();

         /**
         * @brief Destructor
         */
         ~Enstrophy();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_ENSTROPHY
