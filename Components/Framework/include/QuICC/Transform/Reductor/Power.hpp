/**
 * @file Power.hpp
 * @brief Reductor transform operator Reductor::Power
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWER_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWER_HPP

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
    * @brief Reductor transform operator Reductor::Power
    */
   class Power: public IRegisterId<Power>
   {
      public:
         /**
          * @brief Constructor
          */
         Power();

         friend class IRegisterId<Power>;

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWER_HPP
