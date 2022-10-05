/**
 * @file P0.hpp
 * @brief Forward transform operator Forward::P0
 */

#ifndef QUICC_TRANSFORM_FORWARD_P0_HPP
#define QUICC_TRANSFORM_FORWARD_P0_HPP

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
    * @brief Forward transform operator Forward::P0
    */
   class P0: public IRegisterId<P0>
   {
      public:
         /**
          * @brief Constructor
          */
         P0();

         friend class IRegisterId<P0>;

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

#endif // QUICC_TRANSFORM_FORWARD_P0_HPP
