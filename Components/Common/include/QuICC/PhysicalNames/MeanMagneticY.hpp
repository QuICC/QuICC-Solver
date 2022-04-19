/**
 * @file MeanMagneticY.hpp
 * @brief Mean magnetic Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANMAGNETICY_HPP
#define QUICC_PHYSICALNAMES_MEANMAGNETICY_HPP

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
    * @brief Mean magnetic Y physical name
    */
   class MeanMagneticY: public IRegisterId<MeanMagneticY>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanMagneticY();

         friend class IRegisterId<MeanMagneticY>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANMAGNETICY_HPP
