/**
 * @file DxMeanTemperature.hpp
 * @brief D_x mean temperature physical name 
 */

#ifndef QUICC_PHYSICALNAMES_DXMEANTEMPERATURE_HPP
#define QUICC_PHYSICALNAMES_DXMEANTEMPERATURE_HPP

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
    * @brief D_x mean temperature physical name
    */
   class DxMeanTemperature: public IRegisterId<DxMeanTemperature>
   {
      public:
         /**
          * @brief Constructor
          */
         DxMeanTemperature();

         friend class IRegisterId<DxMeanTemperature>;
      
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

#endif // QUICC_PHYSICALNAMES_DXMEANTEMPERATURE_HPP
