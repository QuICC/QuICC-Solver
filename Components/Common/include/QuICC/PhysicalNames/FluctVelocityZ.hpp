/**
 * @file FluctVelocityZ.hpp
 * @brief Fluctuating velocity Z physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTVELOCITYZ_HPP
#define QUICC_PHYSICALNAMES_FLUCTVELOCITYZ_HPP

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
    * @brief Fluctuating velocity Z physical name
    */
   class FluctVelocityZ: public IRegisterId<FluctVelocityZ>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctVelocityZ();

         friend class IRegisterId<FluctVelocityZ>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTVELOCITYZ_HPP
