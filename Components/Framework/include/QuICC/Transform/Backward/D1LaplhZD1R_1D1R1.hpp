/**
 * @file D1LaplhZD1R_1D1R1.hpp
 * @brief Backward projection operator D1LaplhZD1R_1D1R1 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R11_HPP
#define QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R11_HPP

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
    * @brief Backward projection operator D1LaplhZD1R_1D1R1
    */
   class D1LaplhZD1R_1D1R1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D1LaplhZD1R_1D1R1();

         /**
          * @brief Destructor
          */
         virtual ~D1LaplhZD1R_1D1R1();

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

#endif // QUICC_TRANSFORM_BACKWARD_D1LAPLHZD1R_1D1R11_HPP
