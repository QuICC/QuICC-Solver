/**
 * @file D1.hpp
 * @brief Backward projection operator D1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1_HPP
#define QUICC_TRANSFORM_BACKWARD_D1_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Backward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Backward projection operator D1
    */
   class D1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D1();

         /**
          * @brief Destructor
          */
         virtual ~D1();

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

#endif // QUICC_TRANSFORM_BACKWARD_D1_HPP
