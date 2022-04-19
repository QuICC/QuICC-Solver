/**
 * @file VelocityX.hpp
 * @brief Velocity X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VELOCITYX_HPP
#define QUICC_PHYSICALNAMES_VELOCITYX_HPP

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
    * @brief Velocity X physical name
    */
   class VelocityX: public IRegisterId<VelocityX>
   {
      public:
         /**
          * @brief Constructor
          */
         VelocityX();

         friend class IRegisterId<VelocityX>;
      
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

#endif // QUICC_PHYSICALNAMES_VELOCITYX_HPP
