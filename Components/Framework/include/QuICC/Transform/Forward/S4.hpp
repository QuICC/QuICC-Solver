/**
 * @file S4.hpp
 * @brief Forward projection operator S4 
 */

#ifndef QUICC_TRANSFORM_FORWARD_S4_HPP
#define QUICC_TRANSFORM_FORWARD_S4_HPP

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
    * @brief Forward projection operator S4
    */
   class S4: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         S4();

         /**
          * @brief Destructor
          */
         virtual ~S4();

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

#endif // QUICC_TRANSFORM_FORWARD_S4_HPP
