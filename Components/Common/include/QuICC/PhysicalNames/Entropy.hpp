/**
 * @file Entropy.hpp
 * @brief Entropy physical name 
 */

#ifndef QUICC_PHYSICALNAMES_ENTROPY_HPP
#define QUICC_PHYSICALNAMES_ENTROPY_HPP

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
    * @brief Entropy physical name
    */
   class Entropy: public IRegisterId<Entropy>
   {
      public:
         /**
          * @brief Constructor
          */
         Entropy();

         friend class IRegisterId<Entropy>;
      
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

#endif // QUICC_PHYSICALNAMES_ENTROPY_HPP
