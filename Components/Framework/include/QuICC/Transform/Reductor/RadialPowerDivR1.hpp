/**
 * @file RadialPowerDivR1.hpp
 * @brief Reductor transform operator Reductor::RadialPowerDivR1
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1_HPP
#define QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1_HPP

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
    * @brief Reductor transform operator Reductor::RadialPowerDivR1
    */
   class RadialPowerDivR1: public IRegisterId<RadialPowerDivR1>
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPowerDivR1();

         friend class IRegisterId<RadialPowerDivR1>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1_HPP
