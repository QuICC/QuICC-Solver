/**
 * @file FluctMagneticZ.hpp
 * @brief Fluctuating magnetic Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTMAGNETICZ_HPP
#define QUICC_PHYSICALNAMES_FLUCTMAGNETICZ_HPP

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
    * @brief Fluctuating magnetic Z physical name
    */
   class FluctMagneticZ: public IRegisterId<FluctMagneticZ>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctMagneticZ();

         friend class IRegisterId<FluctMagneticZ>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTMAGNETICZ_HPP
