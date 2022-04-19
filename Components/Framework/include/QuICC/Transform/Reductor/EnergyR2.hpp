/**
 * @file EnergyR2.hpp
 * @brief Reductor operator EnergyR2 
 */

#ifndef QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP
#define QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP

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
    * @brief Reductor projection operator EnergyR2
    */
   class EnergyR2: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergyR2();

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

#endif // QUICC_TRANSFORM_REDUCTOR_ENERGYR2_HPP
