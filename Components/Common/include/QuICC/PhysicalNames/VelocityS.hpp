/**
 * @file VelocityS.hpp
 * @brief Velocity S physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VELOCITYS_HPP
#define QUICC_PHYSICALNAMES_VELOCITYS_HPP

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
    * @brief Velocity S physical name
    */
   class VelocityS: public IRegisterId<VelocityS>
   {
      public:
         /**
          * @brief Constructor
          */
         VelocityS();

         friend class IRegisterId<VelocityS>;
      
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

#endif // QUICC_PHYSICALNAMES_VELOCITYS_HPP
