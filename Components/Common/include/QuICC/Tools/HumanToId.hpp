/** 
 * @file HumanToId.hpp
 * @brief Small routines to convert human string to enum ID
 */

#ifndef QUICC_TOOLS_HUMANTOID_HPP
#define QUICC_TOOLS_HUMANTOID_HPP

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
    *  @brief Small routines to convert human string to enum ID
    */
   struct HumanToId
   {
      /**
       * @brief Convert string to field component ID
       */
      static FieldComponents::Spectral::Id toComp(const std::string& id);

      private:
         /**
          * @brief Constructor
          */
         HumanToId() = default;

         /**
          * @brief Constructor
          */
         ~HumanToId() = default;
   };
}
}

#endif // QUICC_TOOLS_HUMANTOID_HPP
