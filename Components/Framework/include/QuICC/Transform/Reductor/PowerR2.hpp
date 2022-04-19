/**
 * @file PowerR2.hpp
 * @brief Reductor operator PowerR2 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP

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
    * @brief Reductor projection operator PowerR2
    */
   class PowerR2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         PowerR2();

         /**
          * @brief Destructor
          */
         virtual ~PowerR2();

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERR2_HPP
