/**
 * @file D1ZP0.hpp
 * @brief Forward projection operator D1ZP0 
 */

#ifndef QUICC_TRANSFORM_FORWARD_D1ZP0_HPP
#define QUICC_TRANSFORM_FORWARD_D1ZP0_HPP

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
    * @brief Forward projection operator D1ZP0
    */
   class D1ZP0: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D1ZP0();

         /**
          * @brief Destructor
          */
         virtual ~D1ZP0();

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

#endif // QUICC_TRANSFORM_FORWARD_D1ZP0_HPP
