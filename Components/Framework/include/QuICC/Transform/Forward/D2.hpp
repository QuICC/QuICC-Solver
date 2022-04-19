/**
 * @file D2.hpp
 * @brief Forward projection operator D2 
 */

#ifndef QUICC_TRANSFORM_FORWARD_D2_HPP
#define QUICC_TRANSFORM_FORWARD_D2_HPP

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
    * @brief Forward projection operator D2
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

#endif // QUICC_TRANSFORM_FORWARD_D2_HPP
