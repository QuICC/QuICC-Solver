/**
 * @file VorticityZ.hpp
 * @brief Vorticity Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VORTICITYZ_HPP
#define QUICC_PHYSICALNAMES_VORTICITYZ_HPP

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
    * @brief Vorticity Z physical name
    */
   class VorticityZ: public IRegisterId<VorticityZ>
   {
      public:
         /**
          * @brief Constructor
          */
         VorticityZ();

         friend class IRegisterId<VorticityZ>;
      
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

#endif // QUICC_PHYSICALNAMES_VORTICITYZ_HPP
