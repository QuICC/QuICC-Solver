/**
 * @file MeanVelocityY.hpp
 * @brief Mean velocity Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANVELOCITYY_HPP
#define QUICC_PHYSICALNAMES_MEANVELOCITYY_HPP

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
    * @brief Mean velocity Y physical name
    */
   class MeanVelocityY: public IRegisterId<MeanVelocityY>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanVelocityY();

         friend class IRegisterId<MeanVelocityY>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANVELOCITYY_HPP
