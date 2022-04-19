/** 
 * @file Energy.hpp
 * @brief Definitions and names use by the energy writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_ENERGY
#define QUICC_IO_VARIABLE_TAGS_ENERGY

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
    * @brief Definitions and names use by the energy writers
    */
   class Energy
   {
      public:
         /**
          * @brief HEADER part for energy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of energy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Energy();

         /**
         * @brief Destructor
         */
         ~Energy();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_ENERGY
