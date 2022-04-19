/**
 * @file Sin_1Dphi.hpp
 * @brief Forward projection operator Sin_1Dphi 
 */

#ifndef QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP
#define QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP

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
    * @brief Forward projection operator Sin_1Dphi
    */
   class Sin_1Dphi: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Sin_1Dphi();

         /**
          * @brief Destructor
          */
         virtual ~Sin_1Dphi();

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

#endif // QUICC_TRANSFORM_FORWARD_SIN_1DPHI_HPP
