/**
 * @file R_1Pm.hpp
 * @brief Backward projection operator R_1Pm 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1PM_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1PM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Backward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Backward projection operator R_1Pm
    */
   class R_1Pm: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         R_1Pm();

         /**
          * @brief Destructor
          */
         virtual ~R_1Pm();

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1PM_HPP
