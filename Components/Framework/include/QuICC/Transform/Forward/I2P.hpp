/**
 * @file I2P.hpp
 * @brief Forward transform operator Forward::I2P
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2P_HPP
#define QUICC_TRANSFORM_FORWARD_I2P_HPP

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
    * @brief Forward transform operator Forward::I2P
    */
   class I2P: public IRegisterId<I2P>
   {
      public:
         /**
          * @brief Constructor
          */
         I2P();

         friend class IRegisterId<I2P>;

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

#endif // QUICC_TRANSFORM_FORWARD_I2P_HPP
