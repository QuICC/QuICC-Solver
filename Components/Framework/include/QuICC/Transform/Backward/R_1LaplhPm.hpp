/**
 * @file R_1LaplhPm.hpp
 * @brief Backward projection operator R_1LaplhPm 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP
#define QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP

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
    * @brief Backward projection operator R_1LaplhPm
    */
   class R_1LaplhPm: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         R_1LaplhPm();

         /**
          * @brief Destructor
          */
         virtual ~R_1LaplhPm();

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

#endif // QUICC_TRANSFORM_BACKWARD_R_1LAPLHPM_HPP
