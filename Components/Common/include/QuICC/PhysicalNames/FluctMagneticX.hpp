/**
 * @file FluctMagneticX.hpp
 * @brief Fluctuating magnetic X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTMAGNETICX_HPP
#define QUICC_PHYSICALNAMES_FLUCTMAGNETICX_HPP

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
    * @brief Fluctuating magnetic X physical name
    */
   class FluctMagneticX: public IRegisterId<FluctMagneticX>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctMagneticX();

         friend class IRegisterId<FluctMagneticX>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTMAGNETICX_HPP
