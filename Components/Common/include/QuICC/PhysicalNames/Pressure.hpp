/**
 * @file Pressure.hpp
 * @brief Pressure physical name 
 */

#ifndef QUICC_PHYSICALNAMES_PRESSURE_HPP
#define QUICC_PHYSICALNAMES_PRESSURE_HPP

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
    * @brief Pressure physical name
    */
   class Pressure: public IRegisterId<Pressure>
   {
      public:
         /**
          * @brief Constructor
          */
         Pressure();

         friend class IRegisterId<Pressure>;
      
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

#endif // QUICC_PHYSICALNAMES_PRESSURE_HPP
