/**
 * @file PowerSLAPLR2.hpp
 * @brief Reductor transform operator Reductor::PowerSLAPLR2
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP

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
    * @brief Reductor transform operator Reductor::PowerSLAPLR2
    */
   class PowerSLAPLR2: public IRegisterId<PowerSLAPLR2>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerSLAPLR2();

         friend class IRegisterId<PowerSLAPLR2>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP
