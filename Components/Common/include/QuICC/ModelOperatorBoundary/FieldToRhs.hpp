/**
 * @file FieldToRhs.hpp
 * @brief FieldToRhs ModelOperatorBoundary
 */

#ifndef QUICC_MODELOPERATORBOUNDARY_FIELDTORHS_HPP
#define QUICC_MODELOPERATORBOUNDARY_FIELDTORHS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperatorBoundary/IRegisterId.hpp"

namespace QuICC {

namespace ModelOperatorBoundary {

   /**
    * @brief FieldToRhs ModelOperatorBoundary
    */
   class FieldToRhs: public IRegisterId<FieldToRhs>
   {
      public:
         /**
          * @brief Constructor
          */
         FieldToRhs();

         friend class IRegisterId<FieldToRhs>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   };

}
}

#endif // QUICC_MODELOPERATORBOUNDARY_FIELDTORHS_HPP
