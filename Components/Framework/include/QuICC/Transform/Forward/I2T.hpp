/**
 * @file I2T.hpp
 * @brief Forward transform operator Forward::I2T
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2T_HPP
#define QUICC_TRANSFORM_FORWARD_I2T_HPP

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
    * @brief Forward transform operator Forward::I2T
    */
   class I2T: public IRegisterId<I2T>
   {
      public:
         /**
          * @brief Constructor
          */
         I2T();

         friend class IRegisterId<I2T>;

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

#endif // QUICC_TRANSFORM_FORWARD_I2T_HPP
