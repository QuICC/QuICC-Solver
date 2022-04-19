/**
 * @file D3.hpp
 * @brief Backward projection operator D3 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D3_HPP
#define QUICC_TRANSFORM_BACKWARD_D3_HPP

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
    * @brief Backward projection operator D3
    */
   class D3: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D3();

         /**
          * @brief Destructor
          */
         virtual ~D3();

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

#endif // QUICC_TRANSFORM_BACKWARD_D3_HPP
