/**
 * @file MeanVelocityZ.hpp
 * @brief Mean velocity Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANVELOCITYZ_HPP
#define QUICC_PHYSICALNAMES_MEANVELOCITYZ_HPP

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
    * @brief Mean velocity Z physical name
    */
   class MeanVelocityZ: public IRegisterId<MeanVelocityZ>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanVelocityZ();

         friend class IRegisterId<MeanVelocityZ>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANVELOCITYZ_HPP
