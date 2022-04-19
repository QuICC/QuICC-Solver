/**
 * @file Temperature.hpp
 * @brief Temperature physical name 
 */

#ifndef QUICC_PHYSICALNAMES_TEMPERATURE_HPP
#define QUICC_PHYSICALNAMES_TEMPERATURE_HPP

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
    * @brief Temperature physical name
    */
   class Temperature: public IRegisterId<Temperature>
   {
      public:
         /**
          * @brief Constructor
          */
         Temperature();

         friend class IRegisterId<Temperature>;
      
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

#endif // QUICC_PHYSICALNAMES_TEMPERATURE_HPP
