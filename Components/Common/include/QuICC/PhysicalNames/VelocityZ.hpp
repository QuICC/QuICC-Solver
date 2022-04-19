/**
 * @file VelocityZ.hpp
 * @brief Velocity Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VELOCITYZ_HPP
#define QUICC_PHYSICALNAMES_VELOCITYZ_HPP

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
    * @brief Velocity Z physical name
    */
   class VelocityZ: public IRegisterId<VelocityZ>
   {
      public:
         /**
          * @brief Constructor
          */
         VelocityZ();

         friend class IRegisterId<VelocityZ>;
      
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

#endif // QUICC_PHYSICALNAMES_VELOCITYZ_HPP
