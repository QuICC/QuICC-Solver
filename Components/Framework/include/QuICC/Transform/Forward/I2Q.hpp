/**
 * @file I2Q.hpp
 * @brief Forward transform operator Forward::I2Q
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2Q_HPP
#define QUICC_TRANSFORM_FORWARD_I2Q_HPP

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
    * @brief Forward transform operator Forward::I2Q
    */
   class I2Q: public IRegisterId<I2Q>
   {
      public:
         /**
          * @brief Constructor
          */
         I2Q();

         friend class IRegisterId<I2Q>;

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

#endif // QUICC_TRANSFORM_FORWARD_I2Q_HPP
