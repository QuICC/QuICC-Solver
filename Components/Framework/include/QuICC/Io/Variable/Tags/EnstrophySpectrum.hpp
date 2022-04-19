/** 
 * @file Enstrophy.hpp
 * @brief Definitions and names use by the enstrophy spectrum writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_ENSTROPHYSPECTRUM
#define QUICC_IO_VARIABLE_TAGS_ENSTROPHYSPECTRUM

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
    * @brief Definitions and names use by the enstrophy spectrum writers
    */
   class EnstrophySpectrum
   {
      public:
         /**
          * @brief HEADER part for enstrophy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for enstrophy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of enstrophy files
          */
         static const std::string   BASENAME;

         /**
          * @brief BASENAME of L spectrum files
          */
         static const std::string   LBASENAME;

         /**
          * @brief BASENAME of M spectrum files
          */
         static const std::string   MBASENAME;

         /**
          * @brief EXTENSION of enstrophy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         EnstrophySpectrum();

         /**
         * @brief Destructor
         */
         ~EnstrophySpectrum();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_ENSTROPHYSPECTRUM
