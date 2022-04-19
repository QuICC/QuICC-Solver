/**
 * @file T.hpp
 * @brief Forward projection operator T 
 */

#ifndef QUICC_TRANSFORM_FORWARD_T_HPP
#define QUICC_TRANSFORM_FORWARD_T_HPP

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
    * @brief Forward projection operator T
    */
   class T: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         T();

         /**
          * @brief Destructor
          */
         virtual ~T();

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

#endif // QUICC_TRANSFORM_FORWARD_T_HPP
