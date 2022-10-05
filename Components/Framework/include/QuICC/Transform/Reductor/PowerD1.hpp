/**
 * @file PowerD1.hpp
 * @brief Reductor transform operator Reductor::PowerD1
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERD1_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERD1_HPP

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
    * @brief Reductor transform operator Reductor::PowerD1
    */
   class PowerD1: public IRegisterId<PowerD1>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1();

         friend class IRegisterId<PowerD1>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERD1_HPP
