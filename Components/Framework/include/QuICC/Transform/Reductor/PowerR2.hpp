/**
 * @file PowerR2.hpp
 * @brief Reductor transform operator Reductor::PowerR2
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP

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
    * @brief Reductor transform operator Reductor::PowerR2
    */
   class PowerR2: public IRegisterId<PowerR2>
   {
      public:
         /**
          * @brief Constructor
          */
         PowerR2();

         friend class IRegisterId<PowerR2>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP
