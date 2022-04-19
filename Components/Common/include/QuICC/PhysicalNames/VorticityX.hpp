/**
 * @file VorticityX.hpp
 * @brief Vorticity X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VORTICITYX_HPP
#define QUICC_PHYSICALNAMES_VORTICITYX_HPP

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
    * @brief Vorticity X physical name
    */
   class VorticityX: public IRegisterId<VorticityX>
   {
      public:
         /**
          * @brief Constructor
          */
         VorticityX();

         friend class IRegisterId<VorticityX>;
      
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

#endif // QUICC_PHYSICALNAMES_VORTICITYX_HPP
