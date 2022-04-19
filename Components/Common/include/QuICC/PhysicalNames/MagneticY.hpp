/**
 * @file MagneticY.hpp
 * @brief Magnetic Y physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MAGNETICY_HPP
#define QUICC_PHYSICALNAMES_MAGNETICY_HPP

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
    * @brief Magnetic Y physical name
    */
   class MagneticY: public IRegisterId<MagneticY>
   {
      public:
         /**
          * @brief Constructor
          */
         MagneticY();

         friend class IRegisterId<MagneticY>;
      
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

#endif // QUICC_PHYSICALNAMES_MAGNETICY_HPP
