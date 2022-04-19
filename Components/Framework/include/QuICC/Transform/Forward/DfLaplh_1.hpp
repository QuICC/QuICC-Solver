/**
 * @file DfLaplh_1.hpp
 * @brief Forward projection operator DfLaplh_1 
 */

#ifndef QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP
#define QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Forward/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   /**
    * @brief Forward projection operator DfLaplh_1
    */
   class DfLaplh_1: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         DfLaplh_1();

         /**
          * @brief Destructor
          */
         virtual ~DfLaplh_1();

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

#endif // QUICC_TRANSFORM_FORWARD_DFLAPLH_1_HPP
