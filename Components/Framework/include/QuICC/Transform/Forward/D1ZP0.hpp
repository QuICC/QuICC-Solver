/**
 * @file D1ZP0.hpp
 * @brief Forward transform operator Forward::D1ZP0
 */

#ifndef QUICC_TRANSFORM_FORWARD_D1ZP0_HPP
#define QUICC_TRANSFORM_FORWARD_D1ZP0_HPP

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
    * @brief Forward transform operator Forward::D1ZP0
    */
   class D1ZP0: public IRegisterId<D1ZP0>
   {
      public:
         /**
          * @brief Constructor
          */
         D1ZP0();

         friend class IRegisterId<D1ZP0>;

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

#endif // QUICC_TRANSFORM_FORWARD_D1ZP0_HPP
