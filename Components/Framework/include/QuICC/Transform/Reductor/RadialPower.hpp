/**
 * @file RadialPower.hpp
 * @brief Reductor operator RadialPower
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP
#define QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP

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
    * @brief Reductor projection operator RadialPower
    */
   class RadialPower: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         RadialPower();

         /**
          * @brief Destructor
          */
         virtual ~RadialPower();

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

#endif // QUICC_TRANSFORM_REDUCTOR_RADIALPOWER_HPP
