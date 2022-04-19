/** 
 * @file Nusselt.hpp
 * @brief Definitions and names use by the Nusselt writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_NUSSELT
#define QUICC_IO_VARIABLE_TAGS_NUSSELT

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
    * @brief Definitions and names use by the Nusselt writer
    */
   class Nusselt
   {
      public:
         /**
          * @brief HEADER part for Nusselt file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Nusselt file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Nusselt file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Nusselt file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Nusselt();

         /**
         * @brief Destructor
         */
         ~Nusselt();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_NUSSELT
