/**
 * @file EnergySLAPLR2.hpp
 * @brief Reductor operator EnergySLAPLR2 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP

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
    * @brief Reductor projection operator EnergySLAPLR2
    */
   class EnergySLAPLR2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         EnergySLAPLR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergySLAPLR2();

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYSLAPLR2_HPP
