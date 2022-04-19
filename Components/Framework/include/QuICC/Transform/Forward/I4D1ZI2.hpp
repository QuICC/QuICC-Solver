/**
 * @file I4D1ZI2.hpp
 * @brief Forward projection operator I4D1ZI2 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4D1ZI2_HPP
#define QUICC_TRANSFORM_FORWARD_I4D1ZI2_HPP

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
    * @brief Forward projection operator I4D1ZI2
    */
   class I4D1ZI2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4D1ZI2();

         /**
          * @brief Destructor
          */
         virtual ~I4D1ZI2();

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

#endif // QUICC_TRANSFORM_FORWARD_I4D1ZI2_HPP
