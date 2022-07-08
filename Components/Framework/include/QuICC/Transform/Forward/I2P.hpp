/**
 * @file I2P.hpp
 * @brief Forward projection operator I2P
 */

#ifndef QUICC_TRANSFORM_FORWARD_I2P_HPP
#define QUICC_TRANSFORM_FORWARD_I2P_HPP

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
    * @brief Forward projection operator I2P
    */
   class I2P: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I2P();

         /**
          * @brief Destructor
          */
         virtual ~I2P();

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

#endif // QUICC_TRANSFORM_FORWARD_I2P_HPP
