/**
 * @file Vorticity.hpp
 * @brief Vorticity physical name 
 */

#ifndef QUICC_PHYSICALNAMES_VORTICITY_HPP
#define QUICC_PHYSICALNAMES_VORTICITY_HPP

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
    * @brief Vorticity physical name
    */
   class Vorticity: public IRegisterId<Vorticity>
   {
      public:
         /**
          * @brief Constructor
          */
         Vorticity();

         friend class IRegisterId<Vorticity>;
      
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

#endif // QUICC_PHYSICALNAMES_VORTICITY_HPP
