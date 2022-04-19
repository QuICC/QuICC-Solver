/**
 * @file R_2.hpp
 * @brief Backward projection operator R_2 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_2_HPP
#define QUICC_TRANSFORM_BACKWARD_R_2_HPP

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
    * @brief Backward projection operator R_2
    */
   class R_2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         R_2();

         /**
          * @brief Destructor
          */
         virtual ~R_2();

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

#endif // QUICC_TRANSFORM_BACKWARD_R_2_HPP
