/**
 * @file D1ZP.hpp
 * @brief Backward projection operator D1ZP 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D1ZP_HPP
#define QUICC_TRANSFORM_BACKWARD_D1ZP_HPP

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
    * @brief Backward projection operator D1ZP
    */
   class D1ZP: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D1ZP();

         /**
          * @brief Destructor
          */
         virtual ~D1ZP();

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

#endif // QUICC_TRANSFORM_BACKWARD_D1ZP_HPP
