/**
 * @file Power.hpp
 * @brief Reductor operator Power
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_POWER_HPP
#define QUICC_TRANSFORM_REDUCTOR_POWER_HPP

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
    * @brief Reductor projection operator Power
    */
   class Power: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Power();

         /**
          * @brief Destructor
          */
         virtual ~Power();

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

#endif // QUICC_TRANSFORM_REDUCTOR_POWER_HPP
