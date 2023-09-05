/** 
 * @file MaxAbsoluteFieldValue.hpp
 * @brief Definitions and names use by the MaxAbsoluteFieldValue writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_MAXABSOLUTEFIELDVALUE
#define QUICC_IO_VARIABLE_TAGS_MAXABSOLUTEFIELDVALUE

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
    * @brief Definitions and names use by the MaxAbsoluteFieldValue writer
    */
   class MaxAbsoluteFieldValue
   {
      public:
         /**
          * @brief HEADER part for MaxAbsoluteFieldValue file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for MaxAbsoluteFieldValue file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of MaxAbsoluteFieldValue file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of MaxAbsoluteFieldValue file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         MaxAbsoluteFieldValue();

         /**
         * @brief Destructor
         */
         ~MaxAbsoluteFieldValue();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_MAXABSOLUTEFIELDVALUE
