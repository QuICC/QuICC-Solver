/** 
 * @file IdToHuman.hpp
 * @brief Small routines to convert enum ID into human strings
 */

#ifndef QUICC_TOOLS_IDTOHUMAN_HPP
#define QUICC_TOOLS_IDTOHUMAN_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Tools {

   /**
    * @brief Simple struct to hold the physical names IDs
    */
   struct IdToHuman
   {
      /**
       * @brief Convert ID to string
       */
      static std::string toString(const FieldComponents::Physical::Id id);

      /**
       * @brief Convert ID to string tag
       */
      static std::string toTag(const FieldComponents::Physical::Id id);

      /**
       * @brief Convert ID to string
       */
      static std::string toString(const FieldComponents::Spectral::Id id);

      /**
       * @brief Convert ID to string tag
       */
      static std::string toTag(const FieldComponents::Spectral::Id id);

      private:
         /**
          * @brief Constructor
          */
         IdToHuman() = default;

         /**
          * @brief Constructor
          */
         ~IdToHuman() = default;
   };
}
}

#endif // QUICC_TOOLS_IDTOHUMAN_HPP
