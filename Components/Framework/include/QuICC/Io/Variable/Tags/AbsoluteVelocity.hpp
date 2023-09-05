/** 
 * @file AbsoluteVelocity.hpp
 * @brief Definitions and names use by the AbsoluteVelocity writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_ABSOLUTEVELOCITY
#define QUICC_IO_VARIABLE_TAGS_ABSOLUTEVELOCITY

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
    * @brief Definitions and names use by the AbsoluteVelocity writer
    */
   class AbsoluteVelocity
   {
      public:
         /**
          * @brief HEADER part for AbsoluteVelocity file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for AbsoluteVelocity file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of AbsoluteVelocity file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of AbsoluteVelocity file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         AbsoluteVelocity();

         /**
         * @brief Destructor
         */
         ~AbsoluteVelocity();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_ABSOLUTEVELOCITY
