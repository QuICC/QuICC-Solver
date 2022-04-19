/**
 * @file ImposedMagnetic.hpp
 * @brief ImposedMagnetic physical name 
 */

#ifndef QUICC_PHYSICALNAMES_IMPOSEDMAGNETIC_HPP
#define QUICC_PHYSICALNAMES_IMPOSEDMAGNETIC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalNames/IRegisterId.hpp"

namespace QuICC {

namespace PhysicalNames {

   /**
    * @brief ImposedMagnetic physical name
    */
   class ImposedMagnetic: public IRegisterId<ImposedMagnetic>
   {
      public:
         /**
          * @brief Constructor
          */
         ImposedMagnetic();

         friend class IRegisterId<ImposedMagnetic>;
      
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

#endif // QUICC_PHYSICALNAMES_IMPOSEDMAGNETIC_HPP
