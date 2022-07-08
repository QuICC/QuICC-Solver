/**
 * @file I2T.hpp
 * @brief Forward projection operator I2T
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2T_HPP
#define QUICC_TRANSFORM_FORWARD_I2T_HPP

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
    * @brief Forward projection operator I2T
    */
   class I2T: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2T();

         /**
          * @brief Destructor
          */
         virtual ~I2T();

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

#endif // QUICC_TRANSFORM_FORWARD_I2T_HPP
