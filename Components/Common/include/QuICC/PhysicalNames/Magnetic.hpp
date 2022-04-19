/**
 * @file Magnetic.hpp
 * @brief Magnetic physical name 
 */

#ifndef QUICC_PHYSICALNAMES_MAGNETIC_HPP
#define QUICC_PHYSICALNAMES_MAGNETIC_HPP

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
    * @brief Magnetic physical name
    */
   class Magnetic: public IRegisterId<Magnetic>
   {
      public:
         /**
          * @brief Constructor
          */
         Magnetic();

         friend class IRegisterId<Magnetic>;
      
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

#endif // QUICC_PHYSICALNAMES_MAGNETIC_HPP
