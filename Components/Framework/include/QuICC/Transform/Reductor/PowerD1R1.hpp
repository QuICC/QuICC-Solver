/**
 * @file PowerD1R1.hpp
 * @brief Reductor transform operator Reductor::PowerD1R1
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP

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
    * @brief Reductor transform operator Reductor::PowerD1R1
    */
   class PowerD1R1: public IRegisterId<PowerD1R1>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1();

         friend class IRegisterId<PowerD1R1>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP
