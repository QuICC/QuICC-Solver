/**
 * @file MeanVelocityX.hpp
 * @brief Mean velocity X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANVELOCITYX_HPP
#define QUICC_PHYSICALNAMES_MEANVELOCITYX_HPP

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
    * @brief Mean velocity X physical name
    */
   class MeanVelocityX: public IRegisterId<MeanVelocityX>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanVelocityX();

         friend class IRegisterId<MeanVelocityX>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANVELOCITYX_HPP
