/**
 * @file I2Q.hpp
 * @brief Forward projection operator I2Q
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2Q_HPP
#define QUICC_TRANSFORM_FORWARD_I2Q_HPP

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
    * @brief Forward projection operator I2Q
    */
   class I2Q: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2Q();

         /**
          * @brief Destructor
          */
         virtual ~I2Q();

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

#endif // QUICC_TRANSFORM_FORWARD_I2Q_HPP
