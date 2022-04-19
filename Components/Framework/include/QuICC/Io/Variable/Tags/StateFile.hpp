/** 
 * @file StateFile.hpp
 * @brief Definitions and names use by the state file readers/writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_STATEFILE
#define QUICC_IO_VARIABLE_TAGS_STATEFILE

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
    * @brief Definitions and names use by the state file readers/writers
    */
   class StateFile
   {
      public:
         /**
          * @brief HEADER part for State file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for State file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of State file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of State file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         StateFile();

         /**
         * @brief Destructor
         */
         ~StateFile();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_STATEFILE
