/**
 * @file I2.hpp
 * @brief Forward projection operator I2 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2_HPP
#define QUICC_TRANSFORM_FORWARD_I2_HPP

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
    * @brief Forward projection operator I2
    */
   class I2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2();

         /**
          * @brief Destructor
          */
         virtual ~I2();

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

#endif // QUICC_TRANSFORM_FORWARD_I2_HPP
