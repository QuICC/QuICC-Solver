/** 
 * @file Dipolarity.hpp
 * @brief Definitions and names use by the Dipolarity writer
 */

#ifndef QUICC_IO_VARIABLE_TAGS_DIPOLARITY_HPP
#define QUICC_IO_VARIABLE_TAGS_DIPOLARITY_HPP

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
    * @brief Definitions and names use by the Dipolarity writer
    */
   class Dipolarity
   {
      public:
         /**
          * @brief HEADER part for Dipolarity file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Dipolarity file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Dipolarity file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Dipolarity file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         Dipolarity() = default;

         /**
         * @brief Destructor
         */
         ~Dipolarity() = default;
   };
} // namespace Tags
} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_TAGS_DIPOLARITY_HPP
