/** 
 * @file Continuity.hpp
 * @brief Definitions and names use by the continuity writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_TAGS_CONTINUITY_HPP
#define QUICC_IO_VARIABLE_TAGS_TAGS_CONTINUITY_HPP

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
    * @brief Definitions and names use by the continuity writer
    */
   class Continuity
   {
      public:
         /**
          * @brief HEADER part for Continuity file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Continuity file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Continuity file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Continuity file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Continuity();

         /**
         * @brief Destructor
         */
         ~Continuity();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_TAGS_CONTINUITY_HPP
