/**
 * @file I2D1.hpp
 * @brief Forward transform operator Forward::I2D1
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2D1_HPP
#define QUICC_TRANSFORM_FORWARD_I2D1_HPP

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
    * @brief Forward transform operator Forward::I2D1
    */
   class I2D1: public IRegisterId<I2D1>
   {
      public:
         /**
          * @brief Constructor
          */
         I2D1();

         friend class IRegisterId<I2D1>;

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

#endif // QUICC_TRANSFORM_FORWARD_I2D1_HPP
