/**
 * @file EnergyD1.hpp
 * @brief Reductor transform operator Reductor::EnergyD1
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP

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
    * @brief Reductor transform operator Reductor::EnergyD1
    */
   class EnergyD1: public IRegisterId<EnergyD1>
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1();

         friend class IRegisterId<EnergyD1>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP
