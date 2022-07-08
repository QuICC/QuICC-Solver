/**
 * @file I4P.hpp
 * @brief Forward projection operator I4P
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4P_HPP
#define QUICC_TRANSFORM_FORWARD_I4P_HPP

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
    * @brief Forward projection operator I4P
    */
   class I4P: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4P();

         /**
          * @brief Destructor
          */
         virtual ~I4P();

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

#endif // QUICC_TRANSFORM_FORWARD_I4P_HPP
