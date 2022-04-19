/**
 * @file PowerD1R1.hpp
 * @brief Reductor operator PowerD1R1 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Reductor/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   /**
    * @brief Reductor projection operator PowerD1R1
    */
   class PowerD1R1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         PowerD1R1();

         /**
          * @brief Destructor
          */
         virtual ~PowerD1R1();

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERD1R1_HPP
