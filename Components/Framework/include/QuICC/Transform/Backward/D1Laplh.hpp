/**
 * @file D1Laplh.hpp
 * @brief Backward projection operator D1Laplh 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP

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
    * @brief Backward projection operator D1Laplh
    */
   class D1Laplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D1Laplh();

         /**
          * @brief Destructor
          */
         virtual ~D1Laplh();

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

#endif // QUICC_TRANSFORM_BACKWARD_D1LAPLH_HPP
