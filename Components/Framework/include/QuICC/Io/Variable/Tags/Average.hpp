/** 
 * @file AverageTags.hpp
 * @brief Definitions and names use by the average writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_TAGS_AVERAGE_HPP
#define QUICC_IO_VARIABLE_TAGS_TAGS_AVERAGE_HPP

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
    * @brief Definitions and names use by the average writers
    */
   class Average
   {
      public:
         /**
          * @brief HEADER part for average files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for average files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of average files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of average files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Average();

         /**
         * @brief Destructor
         */
         ~Average();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_TAGS_AVERAGE_HPP
