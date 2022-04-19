/**
 * @file MagneticX.hpp
 * @brief Magnetic X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MAGNETICX_HPP
#define QUICC_PHYSICALNAMES_MAGNETICX_HPP

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
    * @brief Magnetic X physical name
    */
   class MagneticX: public IRegisterId<MagneticX>
   {
      public:
         /**
          * @brief Constructor
          */
         MagneticX();

         friend class IRegisterId<MagneticX>;
      
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

#endif // QUICC_PHYSICALNAMES_MAGNETICX_HPP
