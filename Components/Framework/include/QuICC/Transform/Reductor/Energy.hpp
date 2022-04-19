/**
 * @file Energy.hpp
 * @brief Reductor operator Energy 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP

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
    * @brief Reductor projection operator Energy
    */
   class Energy: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         Energy();

         /**
          * @brief Destructor
          */
         virtual ~Energy();

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGY_HPP
