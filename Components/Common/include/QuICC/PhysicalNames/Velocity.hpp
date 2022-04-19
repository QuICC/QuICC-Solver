/**
 * @file Velocity.hpp
 * @brief Velocity physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VELOCITY_HPP
#define QUICC_PHYSICALNAMES_VELOCITY_HPP

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
    * @brief Velocity physical name
    */
   class Velocity: public IRegisterId<Velocity>
   {
      public:
         /**
          * @brief Constructor
          */
         Velocity();

         friend class IRegisterId<Velocity>;
      
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

#endif // QUICC_PHYSICALNAMES_VELOCITY_HPP
