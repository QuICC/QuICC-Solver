/**
 * @file FluctVelocityY.hpp
 * @brief Fluctuating velocity Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTVELOCITYY_HPP
#define QUICC_PHYSICALNAMES_FLUCTVELOCITYY_HPP

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
    * @brief Fluctuating velocity Y physical name
    */
   class FluctVelocityY: public IRegisterId<FluctVelocityY>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctVelocityY();

         friend class IRegisterId<FluctVelocityY>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTVELOCITYY_HPP
