/**
 * @file Sin_1.hpp
 * @brief Forward projection operator Sin_1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_SIN_1_HPP
#define QUICC_TRANSFORM_FORWARD_SIN_1_HPP

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
    * @brief Forward projection operator Sin_1
    */
   class Sin_1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1();

         /**
          * @brief Destructor
          */
         virtual ~Sin_1();

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

#endif // QUICC_TRANSFORM_FORWARD_SIN_1_HPP
