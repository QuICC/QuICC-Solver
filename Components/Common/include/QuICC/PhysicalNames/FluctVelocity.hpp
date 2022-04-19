/**
 * @file FluctVelocity.hpp
 * @brief Fluctuating velocity physical name 
 */

#ifndef QUICC_PHYSICALNAMES_FLUCTVELOCITY_HPP
#define QUICC_PHYSICALNAMES_FLUCTVELOCITY_HPP

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
    * @brief Fluctuating velocity physical name
    */
   class FluctVelocity: public IRegisterId<FluctVelocity>
   {
      public:
         /**
          * @brief Constructor
          */
         FluctVelocity();

         friend class IRegisterId<FluctVelocity>;
      
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

#endif // QUICC_PHYSICALNAMES_FLUCTVELOCITY_HPP
