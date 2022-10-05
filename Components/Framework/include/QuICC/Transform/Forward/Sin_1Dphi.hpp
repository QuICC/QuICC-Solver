/**
 * @file Sin_1Dphi.hpp
 * @brief Forward transform operator Forward::Sin_1Dphi
 */

#ifndef QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP
#define QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Forward/IRegisterId.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward transform operator Forward::Sin_1Dphi
    */
   class Sin_1Dphi: public IRegisterId<Sin_1Dphi>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1Dphi();

         friend class IRegisterId<Sin_1Dphi>;

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

} // Forward
} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP
