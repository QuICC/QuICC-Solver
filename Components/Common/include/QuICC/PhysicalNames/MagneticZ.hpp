/**
 * @file MagneticZ.hpp
 * @brief Magnetic Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MAGNETICZ_HPP
#define QUICC_PHYSICALNAMES_MAGNETICZ_HPP

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
    * @brief Magnetic Z physical name
    */
   class MagneticZ: public IRegisterId<MagneticZ>
   {
      public:
         /**
          * @brief Constructor
          */
         MagneticZ();

         friend class IRegisterId<MagneticZ>;
      
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

#endif // QUICC_PHYSICALNAMES_MAGNETICZ_HPP
