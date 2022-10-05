/**
 * @file RadialPowerDivR1D1R1.hpp
 * @brief Reductor transform operator Reductor::RadialPowerDivR1D1R1
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP
#define QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP

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
    * @brief Reductor transform operator Reductor::RadialPowerDivR1D1R1
    */
   class RadialPowerDivR1D1R1: public IRegisterId<RadialPowerDivR1D1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPowerDivR1D1R1();

         friend class IRegisterId<RadialPowerDivR1D1R1>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_RADIALPOWERDIVR1D1R1_HPP
