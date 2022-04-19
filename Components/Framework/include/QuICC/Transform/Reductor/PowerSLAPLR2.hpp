/**
 * @file PowerSLAPLR2.hpp
 * @brief Reductor operator PowerSLAPLR2 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP

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
    * @brief Reductor projection operator PowerSLAPLR2
    */
   class PowerSLAPLR2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         PowerSLAPLR2();

         /**
          * @brief Destructor
          */
         virtual ~PowerSLAPLR2();

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWERSLAPLR2_HPP
