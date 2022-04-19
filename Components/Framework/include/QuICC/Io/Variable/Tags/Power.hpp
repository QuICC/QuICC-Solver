/** 
 * @file Power.hpp
 * @brief Definitions and names use by the power writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_POWER
#define QUICC_IO_VARIABLE_TAGS_POWER

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
    * @brief Definitions and names use by the power writers
    */
   class Power
   {
      public:
         /**
          * @brief HEADER part for power files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for power files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of power files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of power files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Power();

         /**
         * @brief Destructor
         */
         ~Power();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_POWER
