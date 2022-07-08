/**
 * @file I2S.hpp
 * @brief Forward projection operator I2S
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2S_HPP
#define QUICC_TRANSFORM_FORWARD_I2S_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Forward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward projection operator I2S
    */
   class I2S: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2S();

         /**
          * @brief Destructor
          */
         virtual ~I2S();

         /**
          * @brief Unique id
          */
         static const std::size_t& id();

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();
   };

}
}
}

#endif // QUICC_TRANSFORM_FORWARD_I2S_HPP
