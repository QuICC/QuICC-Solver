/**
 * @file VorticityY.hpp
 * @brief Vorticity Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VORTICITYY_HPP
#define QUICC_PHYSICALNAMES_VORTICITYY_HPP

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
    * @brief Vorticity Y physical name
    */
   class VorticityY: public IRegisterId<VorticityY>
   {
      public:
         /**
          * @brief Constructor
          */
         VorticityY();

         friend class IRegisterId<VorticityY>;
      
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

#endif // QUICC_PHYSICALNAMES_VORTICITYY_HPP
