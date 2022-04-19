/**
 * @file DfLaplh.hpp
 * @brief Backward projection operator DfLaplh 
 */

#ifndef QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP
#define QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP

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
    * @brief Backward projection operator DfLaplh
    */
   class DfLaplh: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         DfLaplh();

         /**
          * @brief Destructor
          */
         virtual ~DfLaplh();

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

#endif // QUICC_TRANSFORM_BACKWARD_DFLAPLH_HPP
