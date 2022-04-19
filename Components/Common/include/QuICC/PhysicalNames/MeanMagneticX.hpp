/**
 * @file MeanMagneticX.hpp
 * @brief Mean magnetic X physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANMAGNETICX_HPP
#define QUICC_PHYSICALNAMES_MEANMAGNETICX_HPP

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
    * @brief Mean magnetic X physical name
    */
   class MeanMagneticX: public IRegisterId<MeanMagneticX>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanMagneticX();

         friend class IRegisterId<MeanMagneticX>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANMAGNETICX_HPP
