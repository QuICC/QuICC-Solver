/**
 * @file Q.hpp
 * @brief Forward projection operator Q 
 */

#ifndef QUICC_TRANSFORM_FORWARD_Q_HPP
#define QUICC_TRANSFORM_FORWARD_Q_HPP

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
    * @brief Forward projection operator Q
    */
   class Q: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Q();

         /**
          * @brief Destructor
          */
         virtual ~Q();

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

#endif // QUICC_TRANSFORM_FORWARD_Q_HPP
