/**
 * @file FluctTemperature.hpp
 * @brief Fluctuating temperature physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTTEMPERATURE_HPP
#define QUICC_PHYSICALNAMES_FLUCTTEMPERATURE_HPP

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
    * @brief Fluctuating temperature physical name
    */
   class FluctTemperature: public IRegisterId<FluctTemperature>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctTemperature();

         friend class IRegisterId<FluctTemperature>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTTEMPERATURE_HPP
