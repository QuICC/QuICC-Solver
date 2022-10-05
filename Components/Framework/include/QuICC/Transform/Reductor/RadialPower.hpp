/**
 * @file RadialPower.hpp
 * @brief Reductor transform operator Reductor::RadialPower
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP
#define QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP

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
    * @brief Reductor transform operator Reductor::RadialPower
    */
   class RadialPower: public IRegisterId<RadialPower>
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPower();

         friend class IRegisterId<RadialPower>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP
