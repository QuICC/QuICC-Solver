/**
 * @file Sin_1D1Sin.hpp
 * @brief Backward projection operator Sin_1D1Sin 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP
#define QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP

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
    * @brief Backward projection operator Sin_1D1Sin
    */
   class Sin_1D1Sin: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1D1Sin();

         /**
          * @brief Destructor
          */
         virtual ~Sin_1D1Sin();

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

#endif // QUICC_TRANSFORM_BACKWARD_SIN_1D1SIN_HPP
