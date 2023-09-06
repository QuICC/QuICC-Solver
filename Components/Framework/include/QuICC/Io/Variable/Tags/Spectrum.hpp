/** 
 * @file Spectrum.hpp
 * @brief Definitions and names use by the energy spectrum writers
 */

#ifndef QUICC_IO_VARIABLE_TAGS_SPECTRUM
#define QUICC_IO_VARIABLE_TAGS_SPECTRUM

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
    * @brief Definitions and names use by the energy spectrum writers
    */
   class Spectrum
   {
      public:
         /**
          * @brief HEADER part for energy spectrum files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy spectrum files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy spectrum files
          */
         static const std::string   BASENAME;

         /**
          * @brief BASENAME of L energy spectrum files
          */
         static const std::string   LBASENAME;

         /**
          * @brief BASENAME of M energy spectrum files
          */
         static const std::string   MBASENAME;

         /**
          * @brief BASENAME of N energy spectrum files
          */
         static const std::string   NBASENAME;

         /**
          * @brief BASENAME of mode energy spectrum files
          */
         static const std::string   MODEBASENAME;

         /**
          * @brief BASENAME of radial energy spectrum files
          */
         static const std::string   RBASENAME;

         /**
          * @brief EXTENSION of energy spectrum files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Spectrum();

         /**
         * @brief Destructor
         */
         ~Spectrum();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_SPECTRUM
