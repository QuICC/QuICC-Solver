/**
 * @file D2.hpp
 * @brief Backward projection operator D2 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_D2_HPP
#define QUICC_TRANSFORM_BACKWARD_D2_HPP

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
    * @brief Backward projection operator D2
    */
   class D2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         D2();

         /**
          * @brief Destructor
          */
         virtual ~D2();

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

#endif // QUICC_TRANSFORM_BACKWARD_D2_HPP
