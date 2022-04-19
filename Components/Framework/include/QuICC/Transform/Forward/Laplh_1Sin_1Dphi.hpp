/**
 * @file Laplh_1Sin_1Dphi.hpp
 * @brief Forward projection operator Laplh_1Sin_1Dphi 
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1DPHI_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1DPHI_HPP

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
    * @brief Forward projection operator Laplh_1Sin_1Dphi
    */
   class Laplh_1Sin_1Dphi: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Laplh_1Sin_1Dphi();

         /**
          * @brief Destructor
          */
         virtual ~Laplh_1Sin_1Dphi();

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLH_1SIN_1DPHI_HPP
