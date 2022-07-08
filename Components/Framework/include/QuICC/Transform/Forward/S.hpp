/**
 * @file S.hpp
 * @brief Forward projection operator S 
 */

#ifndef QUICC_TRANSFORM_FORWARD_S_HPP
#define QUICC_TRANSFORM_FORWARD_S_HPP

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
    * @brief Forward projection operator S
    */
   class S: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         S();

         /**
          * @brief Destructor
          */
         virtual ~S();

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

#endif // QUICC_TRANSFORM_FORWARD_S_HPP
