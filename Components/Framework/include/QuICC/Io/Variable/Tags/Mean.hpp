/** 
 * @file Mean.hpp
 * @brief Definitions and names use by the Mean writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_MEAN
#define QUICC_IO_VARIABLE_TAGS_MEAN

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
    * @brief Definitions and names use by the Mean writer
    */
   class Mean
   {
      public:
         /**
          * @brief HEADER part for Mean file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Mean file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Mean file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Mean file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Mean();

         /**
         * @brief Destructor
         */
         ~Mean();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_MEAN
