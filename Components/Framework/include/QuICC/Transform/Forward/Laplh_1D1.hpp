/**
 * @file Laplh_1D1.hpp
 * @brief Forward projection operator Laplh_1D1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP

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
    * @brief Forward projection operator Laplh_1D1
    */
   class Laplh_1D1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh_1D1();

         /**
          * @brief Destructor
          */
         virtual ~Laplh_1D1();

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_1D1_HPP
