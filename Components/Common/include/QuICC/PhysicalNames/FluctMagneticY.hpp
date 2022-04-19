/**
 * @file FluctMagneticY.hpp
 * @brief Fluctuating magnetic Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTMAGNETICY_HPP
#define QUICC_PHYSICALNAMES_FLUCTMAGNETICY_HPP

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
    * @brief Fluctuating magnetic Y physical name
    */
   class FluctMagneticY: public IRegisterId<FluctMagneticY>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctMagneticY();

         friend class IRegisterId<FluctMagneticY>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTMAGNETICY_HPP
