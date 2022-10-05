/**
 * @file Sin_1.hpp
 * @brief Forward transform operator Forward::Sin_1
 */

#ifndef QUICC_TRANSFORM_FORWARD_SIN_1_HPP
#define QUICC_TRANSFORM_FORWARD_SIN_1_HPP

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
    * @brief Forward transform operator Forward::Sin_1
    */
   class Sin_1: public IRegisterId<Sin_1>
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1();

         friend class IRegisterId<Sin_1>;

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

#endif // QUICC_TRANSFORM_FORWARD_SIN_1_HPP
