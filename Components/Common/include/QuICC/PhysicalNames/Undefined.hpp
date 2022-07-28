/**
 * @file Undefined.hpp
 * @brief Undefined physical name 
 */

#ifndef QUICC_PHYSICALNAMES_UNDEFINED_HPP
#define QUICC_PHYSICALNAMES_UNDEFINED_HPP

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
    * @brief Undefined physical name
    */
   class Undefined: public IRegisterId<Undefined>
   {
      public:
         /**
          * @brief Constructor
          */
         Undefined();

         friend class IRegisterId<Undefined>;
      
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

#endif // QUICC_PHYSICALNAMES_UNDEFINED_HPP
