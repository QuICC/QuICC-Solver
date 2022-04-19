/**
 * @file Density.hpp
 * @brief Density physical name 
 */

#ifndef QUICC_PHYSICALNAMES_DENSITY_HPP
#define QUICC_PHYSICALNAMES_DENSITY_HPP

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
    * @brief Density physical name
    */
   class Density: public IRegisterId<Density>
   {
      public:
         /**
          * @brief Constructor
          */
         Density();

         friend class IRegisterId<Density>;
      
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

#endif // QUICC_PHYSICALNAMES_DENSITY_HPP
