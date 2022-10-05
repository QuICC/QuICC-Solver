/**
 * @file EnergyR2.hpp
 * @brief Reductor transform operator Reductor::EnergyR2
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Reductor/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   /**
    * @brief Reductor transform operator Reductor::EnergyR2
    */
   class EnergyR2: public IRegisterId<EnergyR2>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2();

         friend class IRegisterId<EnergyR2>;

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

} // Reductor
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP
