/**
 * @file MeanMagneticZ.hpp
 * @brief Mean magnetic Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANMAGNETICZ_HPP
#define QUICC_PHYSICALNAMES_MEANMAGNETICZ_HPP

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
    * @brief Mean magnetic Z physical name
    */
   class MeanMagneticZ: public IRegisterId<MeanMagneticZ>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanMagneticZ();

         friend class IRegisterId<MeanMagneticZ>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANMAGNETICZ_HPP
