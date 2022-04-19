/**
 * @file I4.hpp
 * @brief Forward projection operator I4 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4_HPP
#define QUICC_TRANSFORM_FORWARD_I4_HPP

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
    * @brief Forward projection operator I4
    */
   class I4: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4();

         /**
          * @brief Destructor
          */
         virtual ~I4();

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

#endif // QUICC_TRANSFORM_FORWARD_I4_HPP
