/**
 * @file FluctVelocityX.hpp
 * @brief Fluctuating velocity X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTVELOCITYX_HPP
#define QUICC_PHYSICALNAMES_FLUCTVELOCITYX_HPP

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
    * @brief Fluctuating velocity X physical name
    */
   class FluctVelocityX: public IRegisterId<FluctVelocityX>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctVelocityX();

         friend class IRegisterId<FluctVelocityX>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTVELOCITYX_HPP
