/**
 * @file LaplhD1.hpp
 * @brief Forward projection operator LaplhD1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP
#define QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP

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
    * @brief Forward projection operator LaplhD1
    */
   class LaplhD1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         LaplhD1();

         /**
          * @brief Destructor
          */
         virtual ~LaplhD1();

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

#endif // QUICC_TRANSFORM_FORWARD_LAPLHD1_HPP
