/**
 * @file Sin_1Laplh.hpp
 * @brief Backward projection operator Sin_1Laplh 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_Sin_1Laplh_HPP
#define QUICC_TRANSFORM_BACKWARD_Sin_1Laplh_HPP

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
    * @brief Backward projection operator Sin_1Laplh
    */
   class Sin_1Laplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1Laplh();

         /**
          * @brief Destructor
          */
         virtual ~Sin_1Laplh();

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

#endif // QUICC_TRANSFORM_BACKWARD_Sin_1Laplh_HPP
