/**
 * @file VorticityS.hpp
 * @brief Vorticity S physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VORTICITYS_HPP
#define QUICC_PHYSICALNAMES_VORTICITYS_HPP

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
    * @brief Vorticity S physical name
    */
   class VorticityS: public IRegisterId<VorticityS>
   {
      public:
         /**
          * @brief Constructor
          */
         VorticityS();

         friend class IRegisterId<VorticityS>;
      
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

#endif // QUICC_PHYSICALNAMES_VORTICITYS_HPP
