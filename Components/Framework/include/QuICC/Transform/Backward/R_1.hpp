/**
 * @file R_1.hpp
 * @brief Backward projection operator R_1
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1_HPP

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
    * @brief Backward projection operator R_1
    */
   class R_1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         R_1();

         /**
          * @brief Destructor
          */
         virtual ~R_1();

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1_HPP
