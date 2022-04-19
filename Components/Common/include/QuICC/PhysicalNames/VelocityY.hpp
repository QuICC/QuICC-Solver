/**
 * @file VelocityY.hpp
 * @brief Velocity Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VELOCITYY_HPP
#define QUICC_PHYSICALNAMES_VELOCITYY_HPP

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
    * @brief Velocity Y physical name
    */
   class VelocityY: public IRegisterId<VelocityY>
   {
      public:
         /**
          * @brief Constructor
          */
         VelocityY();

         friend class IRegisterId<VelocityY>;
      
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

#endif // QUICC_PHYSICALNAMES_VELOCITYY_HPP
