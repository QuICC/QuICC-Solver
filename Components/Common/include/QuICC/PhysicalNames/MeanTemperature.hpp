/**
 * @file MeanTemperature.hpp
 * @brief Mean temperature physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MEANTEMPERATURE_HPP
#define QUICC_PHYSICALNAMES_MEANTEMPERATURE_HPP

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
    * @brief Mean temperature physical name
    */
   class MeanTemperature: public IRegisterId<MeanTemperature>
   {
      public:
         /**
          * @brief Constructor
          */
         MeanTemperature();

         friend class IRegisterId<MeanTemperature>;
      
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

#endif // QUICC_PHYSICALNAMES_MEANTEMPERATURE_HPP
