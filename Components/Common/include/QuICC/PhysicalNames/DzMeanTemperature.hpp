/**
 * @file DzMeanTemperature.hpp
 * @brief D_z mean temperature physical name 
 */

#ifndef QUICC_PHYSICALNAMES_DZMEANTEMPERATURE_HPP
#define QUICC_PHYSICALNAMES_DZMEANTEMPERATURE_HPP

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
    * @brief D_z mean temperature physical name
    */
   class DzMeanTemperature: public IRegisterId<DzMeanTemperature>
   {
      public:
         /**
          * @brief Constructor
          */
         DzMeanTemperature();

         friend class IRegisterId<DzMeanTemperature>;
      
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

#endif // QUICC_PHYSICALNAMES_DZMEANTEMPERATURE_HPP
