/**
 * @file MeanMagnetic.hpp
 * @brief Mean magnetic physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANMAGNETIC_HPP
#define QUICC_PHYSICALNAMES_MEANMAGNETIC_HPP

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
    * @brief Mean magnetic physical name
    */
   class MeanMagnetic: public IRegisterId<MeanMagnetic>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanMagnetic();

         friend class IRegisterId<MeanMagnetic>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANMAGNETIC_HPP
