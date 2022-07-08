/**
 * @file I4Q.hpp
 * @brief Forward projection operator I4Q
 */

#ifndef QUICC_TRANSFORM_FORWARD_I4Q_HPP
#define QUICC_TRANSFORM_FORWARD_I4Q_HPP

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
    * @brief Forward projection operator I4Q
    */
   class I4Q: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4Q();

         /**
          * @brief Destructor
          */
         virtual ~I4Q();

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

#endif // QUICC_TRANSFORM_FORWARD_I4Q_HPP
