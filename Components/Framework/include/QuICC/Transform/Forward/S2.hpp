/**
 * @file S2.hpp
 * @brief Forward projection operator S2 
 */

#ifndef QUICC_TRANSFORM_FORWARD_S2_HPP
#define QUICC_TRANSFORM_FORWARD_S2_HPP

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
    * @brief Forward projection operator S2
    */
   class S2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         S2();

         /**
          * @brief Destructor
          */
         virtual ~S2();

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

#endif // QUICC_TRANSFORM_FORWARD_S2_HPP
