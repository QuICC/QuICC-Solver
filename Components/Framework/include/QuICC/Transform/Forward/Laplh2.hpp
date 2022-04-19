/**
 * @file Laplh2.hpp
 * @brief Forward projection operator Laplh2 
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH2_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH2_HPP

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
    * @brief Forward projection operator Laplh2
    */
   class Laplh2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh2();

         /**
          * @brief Destructor
          */
         virtual ~Laplh2();

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH2_HPP
