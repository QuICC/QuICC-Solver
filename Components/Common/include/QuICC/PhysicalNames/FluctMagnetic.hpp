/**
 * @file FluctMagnetic.hpp
 * @brief Fluctuating magnetic physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTMAGNETIC_HPP
#define QUICC_PHYSICALNAMES_FLUCTMAGNETIC_HPP

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
    * @brief Fluctuating magnetic physical name
    */
   class FluctMagnetic: public IRegisterId<FluctMagnetic>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctMagnetic();

         friend class IRegisterId<FluctMagnetic>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTMAGNETIC_HPP
