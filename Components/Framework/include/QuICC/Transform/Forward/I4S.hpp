/**
 * @file I4S.hpp
 * @brief Forward projection operator I4S
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4S_HPP
#define QUICC_TRANSFORM_FORWARD_I4S_HPP

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
    * @brief Forward projection operator I4S
    */
   class I4S: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4S();

         /**
          * @brief Destructor
          */
         virtual ~I4S();

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

#endif // QUICC_TRANSFORM_FORWARD_I4S_HPP
