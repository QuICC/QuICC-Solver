/**
 * @file I6Laplh.hpp
 * @brief Forward projection operator I6Laplh 
 */

#ifndef QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP
#define QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP

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
    * @brief Forward projection operator I6Laplh
    */
   class I6Laplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I6Laplh();

         /**
          * @brief Destructor
          */
         virtual ~I6Laplh();

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

#endif // QUICC_TRANSFORM_FORWARD_I6LAPLH_HPP
