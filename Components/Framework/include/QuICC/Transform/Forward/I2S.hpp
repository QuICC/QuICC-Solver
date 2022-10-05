/**
 * @file I2S.hpp
 * @brief Forward transform operator Forward::I2S
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2S_HPP
#define QUICC_TRANSFORM_FORWARD_I2S_HPP

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
    * @brief Forward transform operator Forward::I2S
    */
   class I2S: public IRegisterId<I2S>
   {
      public:
         /**
          * @brief Constructor
          */
         I2S();

         friend class IRegisterId<I2S>;

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

#endif // QUICC_TRANSFORM_FORWARD_I2S_HPP
