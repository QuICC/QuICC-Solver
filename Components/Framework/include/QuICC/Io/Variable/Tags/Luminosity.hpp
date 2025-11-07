/** 
 * @file Luminosity.hpp
 * @brief Definitions and names use by the Luminosity writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_LUMINOSITY
#define QUICC_IO_VARIABLE_TAGS_LUMINOSITY

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
    * @brief Definitions and names use by the Luminosity writer
    */
   class Luminosity
   {
      public:
         /**
          * @brief HEADER part for Luminosity file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Luminosity file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Luminosity file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Luminosity file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Luminosity();

         /**
         * @brief Destructor
         */
         ~Luminosity();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_LUMINOSITY
