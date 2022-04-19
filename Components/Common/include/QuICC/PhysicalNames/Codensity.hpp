/**
 * @file Codensity.hpp
 * @brief Codensity physical name 
 */

#ifndef QUICC_PHYSICALNAMES_CODENSITY_HPP
#define QUICC_PHYSICALNAMES_CODENSITY_HPP

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
    * @brief Codensity physical name
    */
   class Codensity: public IRegisterId<Codensity>
   {
      public:
         /**
          * @brief Constructor
          */
         Codensity();

         friend class IRegisterId<Codensity>;
      
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

#endif // QUICC_PHYSICALNAMES_CODENSITY_HPP
