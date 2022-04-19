/**
 * @file I2D1.hpp
 * @brief Forward projection operator I2D1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2D1_HPP
#define QUICC_TRANSFORM_FORWARD_I2D1_HPP

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
    * @brief Forward projection operator I2D1
    */
   class I2D1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2D1();

         /**
          * @brief Destructor
          */
         virtual ~I2D1();

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

#endif // QUICC_TRANSFORM_FORWARD_I2D1_HPP
