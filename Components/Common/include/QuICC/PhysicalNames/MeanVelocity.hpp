/**
 * @file MeanVelocity.hpp
 * @brief Mean velocity physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANVELOCITY_HPP
#define QUICC_PHYSICALNAMES_MEANVELOCITY_HPP

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
    * @brief Mean velocity physical name
    */
   class MeanVelocity: public IRegisterId<MeanVelocity>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanVelocity();

         friend class IRegisterId<MeanVelocity>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANVELOCITY_HPP
