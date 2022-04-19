/**
 * @file EnergyD1.hpp
 * @brief Reductor operator EnergyD1 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP

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
    * @brief Reductor projection operator EnergyD1
    */
   class EnergyD1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1();

         /**
          * @brief Destructor
          */
         virtual ~EnergyD1();

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYD1_HPP
