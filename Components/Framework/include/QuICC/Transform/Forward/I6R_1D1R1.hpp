/**
 * @file I6R_1D1R1.hpp
 * @brief Forward projection operator I6R_1D1R1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP
#define QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP

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
    * @brief Forward projection operator I6R_1D1R1
    */
   class I6R_1D1R1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I6R_1D1R1();

         /**
          * @brief Destructor
          */
         virtual ~I6R_1D1R1();

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

#endif // QUICC_TRANSFORM_FORWARD_I6R_1D1R1_HPP
