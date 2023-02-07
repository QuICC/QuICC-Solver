/** 
 * @file FieldProbe.hpp
 * @brief Definitions and names use by the FieldProbe writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_FIELDPROBE
#define QUICC_IO_VARIABLE_TAGS_FIELDPROBE

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
    * @brief Definitions and names use by the FieldProbe writer
    */
   class FieldProbe
   {
      public:
         /**
          * @brief HEADER part for FieldProbe file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for FieldProbe file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of FieldProbe file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of FieldProbe file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         FieldProbe();

         /**
         * @brief Destructor
         */
         ~FieldProbe();
   };

} //Tags
} // Io
} // Variable
} // QuICC

#endif // QUICC_IO_VARIABLE_TAGS_FIELDPROBE
