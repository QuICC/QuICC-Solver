/**
 * @file P0.hpp
 * @brief Backward projection operator P0 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_P0_HPP
#define QUICC_TRANSFORM_BACKWARD_P0_HPP

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
    * @brief Backward projection operator P0
    */
   class P0: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         P0();

         /**
          * @brief Destructor
          */
         virtual ~P0();

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

#endif // QUICC_TRANSFORM_BACKWARD_P0_HPP
