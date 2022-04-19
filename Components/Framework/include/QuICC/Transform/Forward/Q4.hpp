/**
 * @file Q4.hpp
 * @brief Forward projection operator Q4 
 */

#ifndef QUICC_TRANSFORM_FORWARD_Q4_HPP
#define QUICC_TRANSFORM_FORWARD_Q4_HPP

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
    * @brief Forward projection operator Q4
    */
   class Q4: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Q4();

         /**
          * @brief Destructor
          */
         virtual ~Q4();

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

#endif // QUICC_TRANSFORM_FORWARD_Q4_HPP
