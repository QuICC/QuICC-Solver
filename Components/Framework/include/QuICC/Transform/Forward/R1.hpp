/**
 * @file R1.hpp
 * @brief Forward transform operator Forward::R1
 */

#ifndef QUICC_TRANSFORM_FORWARD_R1_HPP
#define QUICC_TRANSFORM_FORWARD_R1_HPP

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
    * @brief Forward transform operator Forward::R1
    */
   class R1: public IRegisterId<R1>
   {
      public:
         /**
          * @brief Constructor
          */
         R1();

         friend class IRegisterId<R1>;

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

#endif // QUICC_TRANSFORM_FORWARD_R1_HPP
